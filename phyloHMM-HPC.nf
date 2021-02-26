nextflow.enable.dsl = 2


params.genomeSize = '1000'
params.recom_len = '200'
params.recom_rate = '0.05'
params.tMRCA = '0.01'
params.nu_sim = '0.2'
params.xml_file = "${PWD}/bin/GTR_template.xml"
params.out =  "${PWD}/results" 


genome = Channel.from(10)
frequencies = Channel.of(' 0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.of('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
//nu_hmm = Channel.from(0.005,0.01,0.02,0.03,0.04)
//mix_prob = Channel.from(0.2,0.3,0.4,0.5,0.6,0.7)
nu_hmm = Channel.from(0.02)
mix_prob = Channel.from(0.7)
repeat_range = Channel.from(1..4)






process BaciSim {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_$filename" }
     maxForks 1
     
     conda "numpy bioconda::dendropy=4.5.2 conda-forge::matplotlib=3.3.4 pandas"
     
    
     input:
         each genome 
         each repeat_range
         each nu_hmm

     output:
         tuple val(repeat_range), path('BaciSimTrees.tree') , emit: BaciSimtrees
         tuple val(repeat_range), path('clonaltree.tree'), emit: clonaltree
         tuple val(repeat_range), path('BaciSim_Log.txt') , emit: recomlog
         tuple val(repeat_range), path('BaciSim_Recombination.jpeg'), emit: SimFig
          
     
     """
       python3.8  $PWD/bin/BaciSim.py -n ${genome} -g ${params.genomeSize} -l ${params.recom_len} -r ${params.recom_rate}  -nu ${params.nu_sim}  -t ${params.tMRCA}
        
     """
}


process seq_gen {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_$filename" }
    maxForks 1

    conda "bioconda::seq-gen"

    input:
        tuple val(repeat_range),  path('BaciSimTrees.tree')
        val f 
        val r 
        
    
    output:
        path "wholegenome.fasta" , emit: wholegenome
    
    """     
     numTrees=\$(wc -l < BaciSimTrees.tree | awk '{ print \$1 }')
     seq-gen  -mGTR  -l${params.genomeSize} -r$r -f$f -s0.2 -of BaciSimTrees.tree  -p \$numTrees > wholegenome.fasta

    """
}


process Gubbins {

    publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
    maxForks 1

    conda  "bioconda::gubbins=2.4.1"


     input:
        each wholegenome
        each repeat_range
        each nu_hmm
        
    
    output:
        path "gubbins.final_tree.tre" , emit: gubbinstree
    
    """   
     run_gubbins.py  -r GTRGAMMA -p gubbins   ${wholegenome} 
      
    """
}


process get_raxml_tree {

    publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }

    conda "bioconda::raxml=8.2.12"


    input:
        path wholegenome
        each repeat_range
        each nu_hmm


    
    output:
        path 'RAxML_bestTree.tree', emit: myRaxML
    
    
    """
     raxmlHPC -m GTRGAMMA   -p 12345 -s ${wholegenome} -N 10 -n tree 
    """
}



process CFML {

   publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }

   conda "bioconda::clonalframeml=1.12"
    
    input:
        path wholegenome
        path myRaxML 
        each repeat_range
        each nu_hmm
        
    
    output:
        path "CFML.labelled_tree.newick" , emit: CFMLtree
        path "CFML.importation_status.txt" , emit: CFML_recom
    
    """ 
    
     ClonalFrameML ${myRaxML} ${wholegenome} CFML >  CFML.result.txt  
      
    """
}


process phyloHMM {

     publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
     
     conda "numpy bioconda::dendropy=4.5.2 conda-forge::matplotlib=3.3.4 pandas scikit-learn conda-forge::hmmlearn"

     input:
        path wholegenome
        path myRaxML 
        path recomlog
        path CFML_recom
        path CFMLtree
        each nu_hmm
//        each mix_prob
        each repeat_range

     output:
        path 'RecomPartial.xml' , emit: partial_XML 
        path 'originalSeq.xml' , emit: original_XML 
        path 'rmse.rmse' , emit : rmse 
        path 'PhyloHMM_Recombination.jpeg' , emit: phyloHMMFig
        path 'CFML_Recombination.jpeg' , emit: CFMLFig

//      -p ${mix_prob} -nu ${nu_hmm}
     """
       python3.8  $PWD/bin/phyloHMM.py -t ${myRaxML} -a ${wholegenome}  -x ${params.xml_file} -l ${recomlog}  -ct ${CFMLtree} -c ${CFML_recom} -nu ${nu_hmm} 
       
        
     """
}



process Beast_alignment {

     publishDir "${params.out}/n_${genome}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_nu_${nu_hmm}_$filename" }
     
     conda "bioconda::beast2"

     input:
         path original_XML
         each genome
         each nu_hmm
         
         
     output:    

         path 'wholegenome.trees' , emit:  beastSeqTree
     
     """
       beast2  ${original_XML}        
     """
}


process treeannotator_alignment {

     publishDir "${params.out}/n_${genome}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_nu_${nu_hmm}_$filename" }
     
     conda "bioconda::beast"

     input:
     
         path beastSeqTree
         each genome
         each nu_hmm        
         
     output:     
         path 'beastSeqTree.nexus' , emit: SeqTree
     
     """
       treeannotator -b 10  ${beastSeqTree}  beastSeqTree.nexus      
     """
}


process convertor_SeqTree {

     publishDir "${params.out}/n_${genome}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_nu_${nu_hmm}_$filename" }
     
     conda "bioconda::dendropy=4.5.2"
     
     input:
        path SeqTree
        each genome
        each nu_hmm
    
     output:
         path 'beastSeqTree.newick' , emit : beastTree
         
     """
       python3.8 $PWD/bin/NexusToNewick.py -t ${SeqTree} -o 'beastSeqTree.newick'
        
     """
}



process Beast_partial {

     publishDir "${params.out}/n_${genome}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_nu_${nu_hmm}_$filename" }

     input:
        path partial_XML
        each genome
        each nu_hmm
         
         
     output:    
         path 'wholegenome.trees' , emit: beastPartialTree
     
     """
       beast  ${partial_XML}        
     """
}




process treeannotator_partial {

     publishDir "${params.out}/n_${genome}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_nu_${nu_hmm}_$filename" }

     input:   
        path beastPartialTree
        each genome
        each nu_hmm
         
         
     output:
     
         path 'beastOurTree.nexus' , emit: beastOurTree

     
     """
       treeannotator -b 10  ${beastPartialTree}  beastOurTree.nexus      
     """
}



process convertor_ourTree {

     publishDir "${params.out}/n_${genome}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_nu_${nu_hmm}_$filename" }

     input:
     
        path beastOurTree
        each genome
        each nu_hmm

    
     output:
         path 'beastHMMTree.newick' , emit:  beastHMMTree
         
     """
       python3.8 $PWD/bin/NexusToNewick.py -t ${beastOurTree} -o 'beastHMMTree.newick'
        
     """
}


process mergeTreeFiles {

    publishDir "${params.out}/n_${genome}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_nu_${nu_hmm}_$filename" }

    input:
         path beastHMMTree
         path myRaxML 
         path beastTree
         path gubbinstree
         path CFMLtree
         each genome
         each nu_hmm
         

    output:
         path 'allOtherTrees.newick' , emit: allOtherTrees
     
     """
       python3.8 $PWD/bin/mergeFiles.py ${beastHMMTree}  ${myRaxML}  ${beastTree}  ${CFMLtree}  ${gubbinstree} > allOtherTrees.newick
       
     """


}


process TreeCmp {

     publishDir "${params.out}/n_${genome}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_nu_${nu_hmm}_$filename" }

     input:
     
         path clonaltree
         path allOtherTrees
         each genome
         each nu_hmm
         
         
     output:
         path 'TreeCmpResult.result' , emit: Comparison     
    
     """
       java -jar treeCmp.jar  -r ${clonaltree}  -i ${allOtherTrees} -d qt pd rf ms um rfw gdu -o TreeCmpResult.result -W
     """
}


process RMSE_summary {


     publishDir "${params.out}/n_${genome}_num_${repeat_range}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_num_${repeat_range}_$filename" }
     

     input: 
         each genome
         path rmse
         each repeat_range

  
     output:
         path 'rmse_summary.csv' , emit: rmse_summary  
         path 'RMSE_comparison.jpeg' , emit: rmse_plot
    
     """
       python3.8 $PWD/bin/rmse_summary.py -f ${rmse}
       
     """
}




workflow {

    BaciSim(genome,repeat_range,nu_hmm)
    
    seq_gen(BaciSim.out.BaciSimtrees,frequencies,rates)
    
//    Gubbins(seq_gen.out.wholegenome,repeat_range,nu_hmm)

//    get_raxml_tree(seq_gen.out.wholegenome,repeat_range,nu_hmm)
//   
//    CFML(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,repeat_range,nu_hmm)
//   
//    phyloHMM(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,BaciSim.out.recomlog,CFML.out.CFML_recom,CFML.out.CFMLtree,nu_hmm,repeat_range)
 
//    collectedFile = phyloHMM.out.rmse.collectFile(name:"collected_rmse.csv",storeDir:"/home/nehleh/work/results/Summary_Results", keepHeader:false) 
  
//    Beast_alignment(phyloHMM.out.original_XML,genome,nu_hmm)
//    
//    treeannotator_alignment(Beast_alignment.out.beastSeqTree,genome,nu_hmm)
//
//    convertor_SeqTree(treeannotator_alignment.out.SeqTree,BaciSim.out.n_genome,nu_hmm)   
//    
//    Beast_partial(phyloHMM.out.partial_XML,BaciSim.out.n_genome,nu_hmm)
//    
//    treeannotator_partial(Beast_partial.out.beastPartialTree,BaciSim.out.n_genome,nu_hmm)
//    
//    convertor_ourTree(treeannotator_partial.out.beastOurTree,BaciSim.out.n_genome,nu_hmm)
//    
//    mergeTreeFiles(convertor_ourTree.out.beastHMMTree,get_raxml_tree.out.myRaxML,convertor_SeqTree.out.beastTree,Gubbins.out.gubbinstree,CFML.out.CFMLtree,BaciSim.out.n_genome,nu_hmm )
//
//    TreeCmp(BaciSim.out.clonaltree,mergeTreeFiles.out.allOtherTrees,BaciSim.out.n_genome,nu_hmm)
//
//    RMSE_summary(BaciSim.out.n_genome,phyloHMM.out.rmse,BaciSim.out.n_repeat)
       
}