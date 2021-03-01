nextflow.enable.dsl = 2


params.genomeSize = '1000'
params.recom_len = '200'
params.recom_rate = '0.05'
params.tMRCA = '0.01'
params.nu_sim = '0.2'
params.xml_file = "${PWD}/bin/GTR_template.xml"
params.out =  "${PWD}/results" 


genome = Channel.value(7)
frequencies = Channel.value(' 0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
nu_hmm = Channel.of(0.02)
mix_prob = Channel.of(0.6)
repeat_range = Channel.value(1)





process BaciSim {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
     maxForks 1
     
     conda "numpy bioconda::dendropy=4.5.2 conda-forge::matplotlib=3.3.4 pandas"
     
    
     input:
         val genome 
         each repeat_range


     output:
         tuple val(repeat_range), path("BaciSimTrees.tree") , emit: BaciSimtrees
         tuple val(repeat_range), path('clonaltree.tree'), emit: clonaltree
         tuple val(repeat_range), path('BaciSim_Log.txt') , emit: recomlog
         tuple val(repeat_range), path('BaciSim_Recombination.jpeg'), emit: SimFig
          
     
     """
       python3.8  $PWD/bin/BaciSim.py -n ${genome} -g ${params.genomeSize} -l ${params.recom_len} -r ${params.recom_rate}  -nu ${params.nu_sim}  -t ${params.tMRCA}
        
     """
}


process seq_gen {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
    maxForks 1   


    conda "bioconda::seq-gen"

    input:

        tuple val(repeat_range),  path('BaciSimTrees.tree')
        val f 
        val r 
          
    output:
        path "wholegenome_${repeat_range}.fasta" , emit: wholegenome
        val repeat_range , emit: range
    
    """ 
     numTrees=\$(wc -l < BaciSimTrees.tree | awk '{ print \$1 }')
     seq-gen  -mGTR  -l${params.genomeSize} -r$r -f$f -s0.2 -of BaciSimTrees.tree -p \$numTrees > wholegenome_${repeat_range}.fasta

    """
}


process Gubbins {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }

    conda  "bioconda::gubbins=2.4.1"


     input:
        path wholegenome
        val repeat_range
        
    
     output:
        path "gubbins.final_tree.tre" , emit: gubbinstree
        
    
     """   
        run_gubbins.py  -r GTRGAMMA -p gubbins   ${wholegenome} 
          
     """
}


process get_raxml_tree {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
    maxForks 1

    conda "bioconda::raxml=8.2.12"

    input:
        path wholegenome
        val repeat_range

    
    output:
        path 'RAxML_bestTree.tree', emit: myRaxML
    
    
    """
     raxmlHPC -m GTRGAMMA   -p 12345 -s ${wholegenome} -N 10 -n tree 
    """
}


process make_xml_seq {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
    conda "numpy bioconda::dendropy"
    

    input:
        path wholegenome
        path myRaxML
        val repeat_range

     output:
        path 'originalSeq.xml' , emit: original_XML 


     """
       python3.8  $PWD/bin/MakeXMLSeq.py -t ${myRaxML} -a ${wholegenome}  -x ${params.xml_file}        
        
     """
}

process Beast_alignment {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
   
     conda "conda-forge::python=3.7  bioconda::beast2"

     input:
         path original_XML
         val repeat_range
         
         
     output:    

         path 'wholegenome.trees' , emit:  beastSeqTree
     
     """
       beast  ${original_XML}        
     """
}


process treeannotator_alignment {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" } 
     
     conda "conda-forge::python=3.7  bioconda::beast2"

     input:
     
         path beastSeqTree 
         val repeat_range     
         
     output:     
         path 'beastSeqTree.nexus' , emit: SeqTree
     
     """
       treeannotator -b 10  ${beastSeqTree}  beastSeqTree.nexus      
     """
}


process convertor_SeqTree {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
     
     conda "conda-forge::python=3.8  bioconda::dendropy"     

     
     input:
        path SeqTree
        val repeat_range
    
     output:
         path 'beastSeqTree.newick' , emit : beastTree
         
     """
       python3.8 $PWD/bin/NexusToNewick.py -t ${SeqTree} -o 'beastSeqTree.newick'
        
     """
}


process CFML {

   publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
   conda "bioconda::clonalframeml=1.12"
    
    input:
        path wholegenome
        path myRaxML 
        val repeat_range
        
    
    output:
        path "CFML.labelled_tree.newick" , emit: CFMLtree
        path "CFML.importation_status.txt" , emit: CFML_recom
    
    """ 
    
     ClonalFrameML ${myRaxML} ${wholegenome} CFML >  CFML.result.txt  
      
    """
}


process CFML_result {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
     conda "numpy bioconda::dendropy=4.5.2 conda-forge::matplotlib=3.3.4 pandas"

     input:
        path wholegenome
        path myRaxML 
        path CFML_recom
        path CFMLtree
        val repeat_range



     output:
        path 'CFML_Recombination.jpeg' , emit: CFMLFig


     """
       python3.8 $PWD/bin/CMFL_result.py -t ${myRaxML} -a ${wholegenome} -ct ${CFMLtree} -c ${CFML_recom} 
       
        
     """
}



process phyloHMM {

     publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
     
     conda "numpy bioconda::dendropy=4.5.2 conda-forge::matplotlib=3.3.4 pandas scikit-learn conda-forge::hmmlearn"

     input:
        path wholegenome
        path myRaxML 
        tuple val(repeat_range), path('recomlog') 
        path CFML_recom
        val repeat_range
        each nu_hmm
        each mix_prob


     output:
        path 'RecomPartial.xml' , emit: partial_XML 
        path 'rmse.csv' , emit : rmse 
        path 'PhyloHMM_Recombination.jpeg' , emit: phyloHMMFig
        path 'Recom_phyloHMM_Log.txt' , emit: phyloHMMLog
        val nu_hmm , emit: my_nu
        val mix_prob , emit: prob
        
        
     """
       python3.8  $PWD/bin/phyloHMM.py -t ${myRaxML} -a ${wholegenome}  -x ${params.xml_file} -l ${recomlog}  -c ${CFML_recom} -nu ${nu_hmm} -p ${mix_prob}
     
     """
}

process RMSE_summary {

     publishDir "$PWD/results/Summary_Results", mode: "copy"
     maxForks 1
     
     conda "conda-forge::matplotlib=3.3.4 pandas"

     input: 
        path collectedRMSE

  
     output:
         path 'RMSE_comparison.jpeg' , emit: rmse_plot
    
     """
       python3.8 $PWD/bin/rmse_summary.py -f rmse.csv
       
     """
}



process Beast_partial {

     publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
     
     conda "conda-forge::python=3.7  bioconda::beast2"



     input:
        path partial_XML
        each genome
        each nu_hmm
        val repeat_range
         
         
     output:    
         path 'wholegenome.trees' , emit: beastPartialTree
     
     """
       beast  ${partial_XML}        
     """
}




process treeannotator_partial {

     publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
     
     conda "conda-forge::python=3.7  bioconda::beast2"

     input:   
        path beastPartialTree
        each genome
        each nu_hmm
        val repeat_range
         
         
     output:
     
         path 'beastOurTree.nexus' , emit: beastOurTree

     
     """
       treeannotator -b 10  ${beastPartialTree}  beastOurTree.nexus      
     """
}



process convertor_ourTree {

     publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
     
     conda "conda-forge::python=3.8 bioconda::dendropy=4.5.2"
     
     input: 
        path beastOurTree
        each genome
        each nu_hmm
        val repeat_range

    
     output:
         path 'beastHMMTree.newick' , emit:  beastHMMTree
         
     """
       python3.8 $PWD/bin/NexusToNewick.py -t ${beastOurTree} -o 'beastHMMTree.newick'
        
     """
}


process mergeTreeFiles {

    publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
    
    maxForks 1

    input:
         path beastHMMTree
         path myRaxML 
         path beastTree
         path gubbinstree
         path CFMLtree
         each genome
         each nu_hmm
         val repeat_range
         

    output:
         path 'allOtherTrees.newick' , emit: allOtherTrees
     
     """
       python3.8 $PWD/bin/mergeFiles.py ${beastHMMTree}  ${myRaxML}  ${beastTree}  ${CFMLtree}  ${gubbinstree} > allOtherTrees.newick
       
     """


}


process TreeCmp {

     publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
     
     maxForks 1

     input:
     
         tuple val(repeat_range), path('clonaltree')
         path allOtherTrees
         each genome
         each nu_hmm
         val repeat_range
         
         
     output:
         path 'TreeCmpResult.result' , emit: Comparison     
    
     """
       java -jar treeCmp.jar  -r ${clonaltree}  -i ${allOtherTrees} -d qt pd rf ms um rfw gdu -o TreeCmpResult.result -W
     """
}


process phyloHMM_beast {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }
     maxForks 1
     
     conda "numpy bioconda::dendropy=4.5.2 conda-forge::matplotlib=3.3.4 pandas scikit-learn conda-forge::hmmlearn"

     input:
        path wholegenome
        path beastHMMTree 
        tuple val(repeat_range), path('recomlog') 
        path CFML_recom
        path CFMLtree
        each nu_hmm
        val repeat_range

//     output:
//        path 'rmse.rmse' , emit : rmse_partialBeast 
//        path 'PhyloHMM_Recombination.jpeg' , emit: phyloHMMFig_partialBeast

     """
       python3.8 $PWD/bin/phyloHMM.py -t ${beastHMMTree} -a ${wholegenome}  -l ${recomlog}  -nu ${nu_hmm} -ct ${CFMLtree} -c ${CFML_recom} -x ${params.xml_file}
       
        
     """
}



workflow {

    BaciSim(genome,repeat_range)
    
    seq_gen(BaciSim.out.BaciSimtrees,frequencies,rates)
    
    Gubbins(seq_gen.out.wholegenome,seq_gen.out.range)

    get_raxml_tree(seq_gen.out.wholegenome,seq_gen.out.range)
    
    make_xml_seq(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,seq_gen.out.range)
    
    Beast_alignment(make_xml_seq.out.original_XML,seq_gen.out.range)

    treeannotator_alignment(Beast_alignment.out.beastSeqTree,seq_gen.out.range)

    convertor_SeqTree(treeannotator_alignment.out.SeqTree,seq_gen.out.range)  
   
    CFML(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,seq_gen.out.range)

    CFML_result(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,CFML.out.CFML_recom,CFML.out.CFMLtree,seq_gen.out.range)
   
    phyloHMM(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,BaciSim.out.recomlog,CFML.out.CFML_recom,seq_gen.out.range,nu_hmm,mix_prob)
    
    collectedRMSE = phyloHMM.out.rmse.collectFile(name:"rmse.csv",storeDir:"${PWD}/results/Summary_Results", keepHeader:false , sort: false) 
    
    RMSE_summary(collectedRMSE)
      
//    Beast_partial(phyloHMM.out.partial_XML,genome,nu_hmm,seq_gen.out.range)
//    
//    treeannotator_partial(Beast_partial.out.beastPartialTree,genome,nu_hmm,seq_gen.out.range)
//    
//    convertor_ourTree(treeannotator_partial.out.beastOurTree,genome,nu_hmm,seq_gen.out.range)
//    
//    mergeTreeFiles(convertor_ourTree.out.beastHMMTree,get_raxml_tree.out.myRaxML,convertor_SeqTree.out.beastTree,Gubbins.out.gubbinstree,CFML.out.CFMLtree,genome,nu_hmm,seq_gen.out.range)
//
//    TreeCmp(BaciSim.out.clonaltree,mergeTreeFiles.out.allOtherTrees,genome,nu_hmm,seq_gen.out.range)
//    
//    phyloHMM_beast(seq_gen.out.wholegenome,convertor_ourTree.out.beastHMMTree,BaciSim.out.recomlog,CFML.out.CFML_recom,CFML.out.CFMLtree,nu_hmm,seq_gen.out.range)

    
       
}