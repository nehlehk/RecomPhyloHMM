nextflow.enable.dsl = 2


params.genomeSize = '1000'
params.recom_len = '500'
params.recom_rate = '0.05'
params.tMRCA = '0.01'
params.nu_sim = '0.2'
params.xml_file = '/home/nehleh/Documents/GTR_template.xml'
params.out = '/home/nehleh/work/results/'


genome = Channel.from(10)
frequencies = Channel.of(' 0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.of('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
//nu_hmm = Channel.from(0.005,0.01,0.02,0.03,0.04)
//mix_prob = Channel.from(0.2,0.3,0.4,0.5,0.6,0.7)
nu_hmm = Channel.from(0.02)
mix_prob = Channel.from(0.7)
repeat_range = Channel.from(1)







process BaciSim {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
     maxForks 1
     
    
     input:
         each genome 
         each repeat_range
         each nu_hmm

     output:
         path "BaciSimTrees.tree", emit: BaciSimtrees
         path "clonaltree.tree", emit: clonaltree
         path "BaciSim_Log.txt" , emit: recomlog
         path "BaciSim_Recombination.jpeg", emit: SimFig
          
     
     """
       python3.8  /home/nehleh/PhyloCode/RecomPhyloHMM/bin/BaciSim.py -n ${genome} -g ${params.genomeSize} -l ${params.recom_len} -r ${params.recom_rate}  -nu ${params.nu_sim}  -t ${params.tMRCA}
        
     """
}


process seq_gen {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
    maxForks 1
    
    

    input:
        path BaciSimtrees 
        val f 
        val r 
        each repeat_range
        each nu_hmm
        
    
    output:
        path "wholegenome.fasta" , emit: wholegenome
    
    """     
     numTrees=\$(wc -l < ${BaciSimtrees} | awk '{ print \$1 }')
     seq-gen  -mGTR  -l${params.genomeSize} -r$r -f$f -s0.2 -of ${BaciSimtrees}  -p \$numTrees > wholegenome.fasta

    """
}


process Gubbins {

    publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
    maxForks 1

     input:
        path wholegenome
        each repeat_range
        each nu_hmm
        
    
    output:
        path "gubbins.final_tree.tre" , emit: gubbinstree
    
    """   
     run_gubbins  -r GTRGAMMA -p gubbins   ${wholegenome} 
      
    """
}


process get_raxml_tree {

    publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }

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
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/phyloHMM.py -t ${myRaxML} -a ${wholegenome}  -x ${params.xml_file} -l ${recomlog}  -ct ${CFMLtree} -c ${CFML_recom} -nu ${nu_hmm} 
       
        
     """
}



process Beast_alignment {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }

     input:
         path original_XML
         each genome
         each nu_hmm
         each repeat_range
         
         
     output:    

         path 'wholegenome.trees' , emit:  beastSeqTree
     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/beast  ${original_XML}        
     """
}


process treeannotator_alignment {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }

     input:
     
         path beastSeqTree
         each genome
         each nu_hmm   
         each repeat_range     
         
     output:     
         path 'beastSeqTree.nexus' , emit: SeqTree
     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/treeannotator -b 10  ${beastSeqTree}  beastSeqTree.nexus      
     """
}


process convertor_SeqTree {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
     
     input:
        path SeqTree
        each genome
        each nu_hmm
        each repeat_range
    
     output:
         path 'beastSeqTree.newick' , emit : beastTree
         
     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/NexusToNewick.py -t ${SeqTree} -o 'beastSeqTree.newick'
        
     """
}



process Beast_partial {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }

     input:
        path partial_XML
        each genome
        each nu_hmm
        each repeat_range
         
         
     output:    
         path 'wholegenome.trees' , emit: beastPartialTree
     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/beast  ${partial_XML}        
     """
}




process treeannotator_partial {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }

     input:   
        path beastPartialTree
        each genome
        each nu_hmm
        each repeat_range
         
         
     output:
     
         path 'beastOurTree.nexus' , emit: beastOurTree

     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/treeannotator -b 10  ${beastPartialTree}  beastOurTree.nexus      
     """
}



process convertor_ourTree {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
     
     
     input:
     
        path beastOurTree
        each genome
        each nu_hmm
        each repeat_range

    
     output:
         path 'beastHMMTree.newick' , emit:  beastHMMTree
         
     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/NexusToNewick.py -t ${beastOurTree} -o 'beastHMMTree.newick'
        
     """
}


process mergeTreeFiles {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }

    input:
         path beastHMMTree
         path myRaxML 
         path beastTree
         path gubbinstree
         path CFMLtree
         each genome
         each nu_hmm
         each repeat_range
         

    output:
         path 'allOtherTrees.newick' , emit: allOtherTrees
     
     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/mergeFiles.py ${beastHMMTree}  ${myRaxML}  ${beastTree}  ${CFMLtree}  ${gubbinstree} > allOtherTrees.newick
       
     """


}


process TreeCmp {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }

     input:
     
         path clonaltree
         path allOtherTrees
         each genome
         each nu_hmm
         each repeat_range
         
         
     output:
         path 'TreeCmpResult.result' , emit: Comparison     
    
     """
       java -jar /home/nehleh/Documents/0_Research/Software/TreeCmp_v2.0-b76/bin/treeCmp.jar  -r ${clonaltree}  -i ${allOtherTrees} -d qt pd rf ms um rfw gdu -o TreeCmpResult.result -W
     """
}


process RMSE_summary {


     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }
     

     input: 
         each genome
         path rmse
         each repeat_range

  
     output:
         path 'rmse_summary.csv' , emit: rmse_summary  
         path 'RMSE_comparison.jpeg' , emit: rmse_plot
    
     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/rmse_summary.py -f ${rmse}
       
     """
}



process phyloHMM_beast {

     publishDir "${params.out}" , mode: 'copy' , saveAs: { filename -> "num_${repeat_range}_nu_${nu_hmm}/num_${repeat_range}_nu_${nu_hmm}_$filename" }

     input:
        path wholegenome
        path beastHMMTree 
        path recomlog
        path CFML_recom
        path CFMLtree
        each nu_hmm
//        each mix_prob
        each repeat_range

     output:
        path 'rmse.rmse' , emit : rmse_partialBeast 
        path 'PhyloHMM_Recombination.jpeg' , emit: phyloHMMFig_partialBeast

//      -p ${mix_prob} -nu ${nu_hmm}
     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/phyloHMM.py -t ${beastHMMTree} -a ${wholegenome}  -l ${recomlog}  -nu ${nu_hmm} -ct ${CFMLtree} -c ${CFML_recom} 
       
        
     """
}





workflow {

    BaciSim(genome,repeat_range,nu_hmm)
    
    seq_gen(BaciSim.out.BaciSimtrees,frequencies,rates,repeat_range,nu_hmm)
    
    Gubbins(seq_gen.out.wholegenome,repeat_range,nu_hmm)

    get_raxml_tree(seq_gen.out.wholegenome,repeat_range,nu_hmm)
   
    CFML(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,repeat_range,nu_hmm)
   
    phyloHMM(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,BaciSim.out.recomlog,CFML.out.CFML_recom,CFML.out.CFMLtree,nu_hmm,repeat_range)
//    nu_hmm,mix_prob,
   
//    collectedFile = phyloHMM.out.rmse.collectFile(name:"collected_rmse.csv",storeDir:"/home/nehleh/work/results/Summary_Results", keepHeader:false) 
  
    Beast_alignment(phyloHMM.out.original_XML,genome,nu_hmm,repeat_range)
    
    treeannotator_alignment(Beast_alignment.out.beastSeqTree,genome,nu_hmm,repeat_range)

    convertor_SeqTree(treeannotator_alignment.out.SeqTree,genome,nu_hmm,repeat_range)   
    
    Beast_partial(phyloHMM.out.partial_XML,genome,nu_hmm,repeat_range)
    
    treeannotator_partial(Beast_partial.out.beastPartialTree,genome,nu_hmm,repeat_range)
    
    convertor_ourTree(treeannotator_partial.out.beastOurTree,genome,nu_hmm,repeat_range)
    
    mergeTreeFiles(convertor_ourTree.out.beastHMMTree,get_raxml_tree.out.myRaxML,convertor_SeqTree.out.beastTree,Gubbins.out.gubbinstree,CFML.out.CFMLtree,genome,nu_hmm,repeat_range)

    TreeCmp(BaciSim.out.clonaltree,mergeTreeFiles.out.allOtherTrees,genome,nu_hmm,repeat_range)
    
    phyloHMM_beast(seq_gen.out.wholegenome,convertor_ourTree.out.beastHMMTree,BaciSim.out.recomlog,CFML.out.CFML_recom,CFML.out.CFMLtree,nu_hmm,repeat_range)

//    RMSE_summary(BaciSim.out.n_genome,phyloHMM.out.rmse,BaciSim.out.n_repeat)
       
}