#!/usr/bin/env nextflow

resultsRoot = params.resultsRoot
etienne_scripts="${workflow.projectDir}/etienne_scripts"

RAW_COUNT_MATRIX = Channel.fromPath( "${params.matrix}" )
CDNA_GTF = Channel.fromPath( "${params.gtf}" )
GTF = Channel.fromPath("${params.gtf}")
SDRF = Channel.fromPath("${params.sdrf}")

matrix_file = file("${params.matrix}")
matrix_name = matrix_file.getSimpleName()

// Read in the .mtx format data

process read_10x {

    conda "${workflow.projectDir}/envs/scanpy.yml"
    
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10

    input:
        file expressionMatrix from RAW_COUNT_MATRIX
        
    output:
        file "${matrix_name}_raw.h5ad" into RAW_ANNDATA

    """
        zipdir=\$(unzip -qql ${expressionMatrix.getBaseName()}.zip | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/||')
        unzip ${expressionMatrix.getBaseName()}
        scanpy-read-10x.py -d \${zipdir}/ -o ${matrix_name}_raw.h5ad -F anndata
    """
}

// Filter cells

process filter_cells {

    conda "${workflow.projectDir}/envs/scanpy.yml"
    
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10

    input:
        file rawData from RAW_ANNDATA
        
    output:
        file "${matrix_name}_filter_cells.h5ad" into FILTER_CELLS_ANNDATA

    """
        scanpy-filter-cells.py -i ${rawData} -p n_genes,n_counts \
            -l ${params.scanpy.filter_cells.min_genes},${params.scanpy.filter_cells.min_counts} \
            -j ${params.scanpy.filter_cells.max_genes},${params.scanpy.filter_cells.max_counts} \
            -o ${matrix_name}_filter_cells.h5ad
    """
}

// Filter genes

process make_genelist {

    conda "${baseDir}/envs/bioconductor-rtracklayer.yml"
    
    memory { 3.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 5

    input:
        file cdnaGtf from CDNA_GTF

    output:
        file 'genes.txt' into GENE_LIST

    """
        #!/usr/bin/env Rscript
        
        suppressPackageStartupMessages(require(rtracklayer))
        annotation <- elementMetadata(import('$cdnaGtf'))
        genes <- unique(annotation[['gene_id']])
        writeLines(genes[ ! is.na(genes)], con = 'genes.txt') 
    """
}


// Extract gene names so we can limit to this list, excluding spikes. We may
// decide to use spikes to normalise in the future, but exclude for now

process filter_genes {

    conda "${workflow.projectDir}/envs/scanpy.yml"
    
    publishDir "$resultsRoot/matrices", mode: 'copy', overwrite: true
    
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10

    input:
        file filterCellsData from FILTER_CELLS_ANNDATA
        file geneList from GENE_LIST        

    output:
        file "${matrix_name}_filter_genes.h5ad" into FILTER_GENES_ANNDATA
        file "${matrix_name}_filter_cells_genes.zip" into FILTER_CELLS_MTX

    """
        mkdir -p ${matrix_name}_filter_cells_genes

        scanpy-filter-genes.py -i ${filterCellsData} -s ${geneList} -p n_cells,n_counts \
            -l ${params.scanpy.filter_genes.min_cells},${params.scanpy.filter_genes.min_counts} \
            -j ${params.scanpy.filter_genes.max_cells},${params.scanpy.filter_genes.max_counts} \
            -o ${matrix_name}_filter_genes.h5ad -x ${matrix_name}_filter_cells_genes/
        
        zip -r ${matrix_name}_filter_cells_genes.zip ${matrix_name}_filter_cells_genes/
        
        rm -rf ${matrix_name}_filter_cells_genes
    """
}

// Normalise expression values across cells

process normalise_data {

    conda "${workflow.projectDir}/envs/scanpy.yml"
    
    publishDir "$resultsRoot/matrices", mode: 'copy', overwrite: true
    
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10

    input:
        file filterGenesData from FILTER_GENES_ANNDATA

    output:
        file "${matrix_name}_normalised.h5ad" into NORMALISED_ANNDATA
        file "${matrix_name}_normalised.zip" into NORMALISED_MTX

    """
        mkdir -p ${matrix_name}_normalised

        scanpy-normalise-data.py -i ${filterGenesData} -s ${params.scanpy.normalise_data.scale_factor} \
             -o ${matrix_name}_normalised.h5ad --save-raw -x ${matrix_name}_normalised/
        
        zip -r ${matrix_name}_normalised.zip ${matrix_name}_normalised/

        rm -rf ${matrix_name}_normalised
    """
}

// Select out the variable genes 

process find_variable_genes {

    conda "${workflow.projectDir}/envs/scanpy.yml"
    
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10

    input:
        file normalisedData from NORMALISED_ANNDATA

    output:
        file "${matrix_name}_findvariablegenes.h5ad" into FIND_VARIABLE_GENES_ANNDATA
        file "variable_genes.png" into VARIABLE_GENES_PLOT

    """
        scanpy-find-variable-genes.py -i ${normalisedData} --flavor ${params.scanpy.find_variable_genes.flavor} \
            -p mean,disp -l ${params.scanpy.find_variable_genes.min_mean},${params.scanpy.find_variable_genes.min_disp} \
            -j ${params.scanpy.find_variable_genes.max_mean},${params.scanpy.find_variable_genes.max_disp} \
            -b ${params.scanpy.find_variable_genes.n_bins} \
            -P variable_genes.png -o ${matrix_name}_findvariablegenes.h5ad

    """
}

// Run scaling (only logarithmize data)

process scale_data {

    conda "${workflow.projectDir}/envs/etienne_scripts.yml"
    
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10

    input:
        file findVariableGenesData from FIND_VARIABLE_GENES_ANNDATA

    output:
        file "${matrix_name}_scaledata.h5ad" into SCALE_DATA_ANNDATA

    """
        python ${etienne_scripts}/scale_data.py ${findVariableGenesData} ${matrix_name}_scaledata.h5ad  
    """
}

// Run principal components analysis

process run_pca {

    conda "${workflow.projectDir}/envs/scanpy.yml"

    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10
    
    publishDir "$resultsRoot/pca", mode: 'copy', overwrite: true
    
    input:
        file scaledData from SCALE_DATA_ANNDATA

    output:
        file "${matrix_name}_pca.h5ad" into PCA_ANNDATA
        file 'embeddings.csv'
        file 'loadings.csv'
        file 'stdev.txt'
        file 'var_ratio.txt'
        file 'pca.png'

    script:

        zero_centre = ''
        if ( params.scanpy.pca.containsKey('zero_centre') && params.scanpy.pca.zero_centre == 'false' ){
            zero_centre = "--no-zero-center"
        }else{
            zero_centre = "--zero-center"
            
        }

        chunked = ''
        if ( params.scanpy.pca.containsKey('chunked') && params.scanpy.pca.chunked == 'true' && params.scanpy.pca.containsKey('chunk_size') ){
            chunked = "--chunked --chunk-size ${params.scanpy.pca.chunk_size}"
        }

        color_by = ''
        if ( params.scanpy.pca.containsKey('color_by') && params.scanpy.pca.color_by != 'none' ){
            color_by = "--color-by ${params.scanpy.pca.color_by}"
        }

        edges = ''
        if ( params.scanpy.pca.containsKey('edges') && params.scanpy.pca.edges != 'false' ){
            edges = '--edges'
        }
        
        use_raw = ''
        if ( params.scanpy.pca.containsKey('use_raw') && params.scanpy.pca.use_raw != 'false' ){
            use_raw = '--use-raw'
        }

        arrows = ''
        if ( params.scanpy.pca.containsKey('arrows') && params.scanpy.pca.arrows != 'false' ){
            arrows = '--arrows'
        }
        
        sort_order = ''
        if ( params.scanpy.pca.containsKey('sort_order') && params.scanpy.pca.sort_order == 'false' ){
            sort_order = '--no-sort-order'
        }
        
        groups = ''
        if ( params.scanpy.pca.containsKey('groups') && params.scanpy.pca.groups != 'none' ){
            groups = "--groups ${params.scanpy.pca.groups}"
        }
  
        frame = ''
        if ( params.scanpy.pca.containsKey('frame') && params.scanpy.pca.frame == 'false' ){
            frame = "--frameoff"
        }
      

        """
            scanpy-run-pca.py -i ${scaledData} -o ${matrix_name}_pca.h5ad \
                --output-embeddings-file embeddings.csv --output-loadings-file loadings.csv \
                --output-stdev-file stdev.txt --output-var-ratio-file var_ratio.txt \
                --n-pcs ${params.scanpy.pca.n_pcs} --svd-solver ${params.scanpy.pca.svd_solver} \
                --random-seed ${params.scanpy.pca.random_seed} ${chunked} --output-plot pca.png \
                $color_by $use_raw $edges $arrows $sort_order --projection ${params.scanpy.pca.projection} \
                --components ${params.scanpy.pca.components} --palette ${params.scanpy.pca.palette} $zero_centre
        """
}

// Filter cells

process neighbours {

    conda "${workflow.projectDir}/envs/scanpy.yml"

    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10
    
    input:
        file pcaData from PCA_ANNDATA
        
    output:
        file "${matrix_name}_neighbours.h5ad" into NEIGHBOURS_ANNDATA

    script:

        knn = ''
        if ( params.scanpy.neighbours.containsKey('knn') && params.scanpy.neighbours.knn == 'false' ){
            knn = "--knn"
        }else{
            knn = "--no-knn"
        }

        use_rep = ''
        if ( params.scanpy.neighbours.containsKey('use_rep') && params.scanpy.neighbours.use_rep != 'none' ){
            use_rep = "--use-rep ${params.scanpy.neighbours.use_rep}"
        }

        """
            scanpy-neighbours.py -i ${pcaData} -o ${matrix_name}_neighbours.h5ad \
                --n-neighbors ${params.scanpy.neighbours.n_neighbours} --n-pcs ${params.scanpy.neighbours.n_pcs} \
                ${use_rep} ${knn} --random-seed ${params.scanpy.neighbours.random_seed} \
                --method ${params.scanpy.neighbours.method} --metric ${params.scanpy.neighbours.metric}
        """
}

// Find clusters

process find_cluster {

    conda "${workflow.projectDir}/envs/scanpy.yml"
    
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10
    
    publishDir "$resultsRoot/clustering", mode: 'copy', overwrite: true

    input:
        file neighboursData from NEIGHBOURS_ANNDATA
        
    output:
        file "${matrix_name}_clusters.h5ad" into CLUSTERS_ANNDATA
        file "clusters.txt"

    script:

        use_weights = ''
        if ( params.scanpy.find_clusters.containsKey('use_weights') && params.scanpy.find_clusters.use_weights != 'false' ){
            use_weights = "--use-weights"
        }

        restrict_to = ''
        if ( params.scanpy.find_clusters.containsKey('restrict_to') && params.scanpy.find_clusters.restrict_to != 'none' ){
            restrict_to = "--restrict-to ${params.scanpy.find_clusters.restrict_to}"
        }

        resolutions = params.scanpy.find_clusters.resolutions.join(",")

        """
            scanpy-find-cluster.py -i ${neighboursData} -o ${matrix_name}_clusters.h5ad \
                --output-text-file clusters.txt --flavor ${params.scanpy.find_clusters.flavor} \
                --resolution ${resolutions} ${restrict_to} ${use_weights} \
                --key-added ${params.scanpy.find_clusters.key_added} --random-seed ${params.scanpy.find_clusters.random_seed}  
        """
}

// Send the clustered object to various downstream operations

CLUSTERS_ANNDATA.into{
    CLUSTERS_ANNDATA_FOR_MARKERS
    CLUSTERS_ANNDATA_FOR_UMAP
    CLUSTERS_ANNDATA_FOR_TSNE
}

// Run UMAP

process run_umap {

    conda "${workflow.projectDir}/envs/scanpy.yml"

    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10
    
    publishDir "$resultsRoot/umap", mode: 'copy', overwrite: true

    input:
        file clusteredAnndata from CLUSTERS_ANNDATA_FOR_UMAP
        
    output:
        file "${matrix_name}_clusters_umap.h5ad" into UMAP_ANNDATA
        file 'embeddings.csv'
        file 'umap.png'

    script:

        maxiter = ''
        if ( params.scanpy.run_umap.containsKey('maxiter') && params.scanpy.run_umap.maxiter != 'none' ){
            maxiter = "--maxiter ${params.scanpy.run_umap.maxiter}"
        }
        
        a = ''
        if ( params.scanpy.run_umap.containsKey('a') && params.scanpy.run_umap.a != 'none' ){
            maxiter = "-a ${params.scanpy.run_umap.a}"
        }
        
        b = ''
        if ( params.scanpy.run_umap.containsKey('b') && params.scanpy.run_umap.b != 'none' ){
            maxiter = "-b ${params.scanpy.run_umap.b}"
        }
    
        use_raw = ''
        if ( params.scanpy.run_umap.containsKey('use_raw') && params.scanpy.run_umap.use_raw != 'false' ){
            use_raw = '--use-raw'
        }

        color_by = ''
        if ( params.scanpy.run_umap.containsKey('color_by') && params.scanpy.run_umap.color_by != 'none' ){
            color_by = "--color-by ${params.scanpy.run_umap.color_by}"
        }

        edges = ''
        if ( params.scanpy.run_umap.containsKey('edges') && params.scanpy.run_umap.edges != 'false' ){
            edges = '--edges'
        }
        
        arrows = ''
        if ( params.scanpy.run_umap.containsKey('arrows') && params.scanpy.run_umap.arrows != 'false' ){
            arrows = '--arrows'
        }
        
        sort_order = ''
        if ( params.scanpy.run_umap.containsKey('sort_order') && params.scanpy.run_umap.sort_order == 'false' ){
            sort_order = '--no-sort-order'
        }
        
        groups = ''
        if ( params.scanpy.run_umap.containsKey('groups') && params.scanpy.run_umap.groups != 'none' ){
            groups = "--groups ${params.scanpy.run_umap.groups}"
        }
        
        """
            scanpy-run-umap.py -i ${clusteredAnndata} -o ${matrix_name}_clusters_umap.h5ad \
                --output-embeddings-file embeddings.csv --min-dist ${params.scanpy.run_umap.min_dist} \
                --spread ${params.scanpy.run_umap.spread} --n-components ${params.scanpy.run_umap.n_components} \
                ${maxiter} --alpha ${params.scanpy.run_umap.alpha} --gamma ${params.scanpy.run_umap.gamma} \
                --negative-sample-rate ${params.scanpy.run_umap.negative_sample_rate} \
                --init-pos ${params.scanpy.run_umap.init_pos} --random-seed ${params.scanpy.run_umap.random_seed} \
                ${a} ${b} --output-plot umap.png ${use_raw} ${color_by} ${edges} ${arrows} ${sort_order} ${groups} \
                --projection ${params.scanpy.run_umap.projection}
        """
}

// Run t-SNE for a range of perplexities

process run_tsne {

    conda "${workflow.projectDir}/envs/scanpy.yml"
    
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10
    
    publishDir "$resultsRoot/tsne", mode: 'copy', overwrite: true

    input:
        file clusteredAnndata from CLUSTERS_ANNDATA_FOR_TSNE
        each perplexity from params.scanpy.run_tsne.perplexities
        
    output:
        file "${matrix_name}_clusters_tsne_${perplexity}.h5ad" into TSNE_ANNDATA
        file "embeddings_${perplexity}.csv"
        file "tsne_${perplexity}.png"

    script:
        
        use_raw = ''
        if ( params.scanpy.run_tsne.containsKey('use_raw') && params.scanpy.run_tsne.use_raw != 'false' ){
            use_raw = '--use-raw'
        }

        color_by = ''
        if ( params.scanpy.run_tsne.containsKey('color_by') && params.scanpy.run_tsne.color_by != 'none' ){
            color_by = "--color-by ${params.scanpy.run_tsne.color_by}"
        }

        edges = ''
        if ( params.scanpy.run_tsne.containsKey('edges') && params.scanpy.run_tsne.edges != 'false' ){
            edges = '--edges'
        }
        
        arrows = ''
        if ( params.scanpy.run_tsne.containsKey('arrows') && params.scanpy.run_tsne.arrows != 'false' ){
            arrows = '--arrows'
        }
        
        sort_order = ''
        if ( params.scanpy.run_tsne.containsKey('sort_order') && params.scanpy.run_tsne.sort_order == 'false' ){
            sort_order = '--no-sort-order'
        }
        
        groups = ''
        if ( params.scanpy.run_tsne.containsKey('groups') && params.scanpy.run_tsne.groups != 'none' ){
            groups = "--groups ${params.scanpy.run_tsne.groups}"
        }
        
        use_rep = ''
        if ( params.scanpy.run_tsne.containsKey('use_rep') && params.scanpy.neighbours.use_rep != 'none' ){
            use_rep = "--use-rep ${params.scanpy.tsne.use_rep}"
        }
        
        frame = ''
        if ( params.scanpy.run_tsne.containsKey('frame') && params.scanpy.run_tsne.frame == 'false' ){
            frame = "--frameoff"
        }

        """
            scanpy-run-tsne.py -i ${clusteredAnndata} -o ${matrix_name}_clusters_tsne_${perplexity}.h5ad \
                --output-embeddings-file embeddings_${perplexity}.csv ${use_rep} --perplexity ${perplexity} \
                --early-exaggeration ${params.scanpy.run_tsne.early_exaggeration} --learning-rate ${params.scanpy.run_tsne.learning_rate} \
                --random-seed ${params.scanpy.run_tsne.random_seed} --output-plot tsne_${perplexity}.png \
                ${use_raw} ${color_by} ${edges} ${arrows} ${sort_order} ${groups} --projection ${params.scanpy.run_tsne.projection} \
                --components ${params.scanpy.run_tsne.components} --palette ${params.scanpy.run_tsne.palette} $frame
        """
}

// Run marker detection 

process find_markers {

    conda "${workflow.projectDir}/envs/scanpy.yml"
    
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "$resultsRoot/markers", mode: 'copy', overwrite: true

    input:
        each resolution from params.scanpy.find_clusters.resolutions
        file clusteredAnndata from CLUSTERS_ANNDATA_FOR_MARKERS
        
    output:
        file "${matrix_name}_clusters_${resolution}_markers.h5ad" into MARKERS_ANNDATA
        file "markers_${resolution}.csv"
        file "markers_${resolution}.png"

    script:
        groups = ''
        if ( params.scanpy.find_markers.containsKey('groups') && params.scanpy.find_markers.groups != 'none' ){
            groups = "--groups ${params.scanpy.find_markers.groups}"
        }

        rankby_abs = ''
        if ( params.scanpy.find_markers.containsKey('rankby_abs') && params.scanpy.find_markers.rankby_abs == 'true' ){
            rankby_abs = "--rankby_abs"
        }
        
        use_raw = ''
        if ( params.scanpy.find_clusters.containsKey('use_raw') && params.scanpy.find_clusters.use_raw == 'false' ){
            use_raw = '--no-raw'
        }

        key = ''
        if ( params.scanpy.find_clusters.containsKey('key') && params.scanpy.find_clusters.key == 'none' ){
            key = "--key ${params.scanpy.find_clusters.key}"
        }
        
        """
            scanpy-find-markers.py -i ${clusteredAnndata} \
                -o "${matrix_name}_clusters_${resolution}_markers.h5ad" \
                --output-text-file markers_${resolution}.csv --groupby "${params.scanpy.find_clusters.key_added}_r${resolution}" ${groups} \
                --reference ${params.scanpy.find_markers.reference} --n-genes ${params.scanpy.find_markers.n_genes} \
                --method ${params.scanpy.find_markers.method} ${rankby_abs} ${use_raw} --output-plot markers_${resolution}.png \
                --show-n-genes ${params.scanpy.find_markers.show_n_genes} ${key} 
        """
}

process merge_matrices {

    conda "${workflow.projectDir}/envs/etienne_scripts.yml"
    publishDir "$resultsRoot",mode: "copy", overwrite: true

    input:
        file sdrf from SDRF
        file gtf from GTF
        file tsne from TSNE_ANNDATA
        file umap from UMAP_ANNDATA

    output:
        file "anndata_merged.h5ad" into RESULTS_SCANPY

    """
        python ${etienne_scripts}/merge_anndata.py -sdrf ${sdrf} -gtf ${gtf} -matrix ${umap} -tsne ${tsne} -o anndata_merged.h5ad  
    """
}
