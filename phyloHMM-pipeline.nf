nextflow.enable.dsl = 2


params.genomeSize = '5000'
params.recom_len = '500'
params.recom_rate = '0.05'
params.nu_sim = '0.2'
params.xml_file = '/home/nehleh/Documents/GTR_template.xml'
params.out = '/home/nehleh/work/results/'


genome = Channel.from(10)
frequencies = Channel.of(' 0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.of('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
nu_hmm = Channel.from(0.02,0.03,0.04,0.05)
mix_prob = Channel.from(0.2,0.3,0.4,0.5,0.6,0.7)
repeat_range = Channel.from( 1..50 )







process BaciSim {

     publishDir "${params.out}/n_${genome}_num_${repeat_range}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_num_${repeat_range}_$filename" }   
     
     input:
         each genome 
         each repeat_range

     output:
         path "BaciSimTrees.tree", emit: BaciSimtrees
         path "clonaltree.tree", emit: clonaltree
         path "BaciSim_Log.txt" , emit: recomlog
         path "BaciSim_Recombination.jpeg", emit: SimFig
         val "${genome}" , emit: n_genome
         val "${repeat_range}" , emit: n_repeat
          
     
     """
       python3.8  /home/nehleh/PhyloCode/RecomPhyloHMM/bin/BaciSim.py -n ${genome} -g ${params.genomeSize} -l ${params.recom_len} -r ${params.recom_rate}  -nu ${params.nu_sim}  
        
     """
}


process seq_gen {

    publishDir "${params.out}/n_${genome}_num_${repeat_range}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_num_${repeat_range}_$filename" }

    input:
        file trees 
        val f 
        val r 
        each genome
        each repeat_range
        
    
    output:
        path "wholegenome.fasta" , emit: wholegenome
    
    """     
     numTrees=\$(wc -l < ${trees} | awk '{ print \$1 }')
     seq-gen  -mGTR  -l${params.genomeSize} -r$r -f$f -s0.2 -of ${trees}  -p \$numTrees > wholegenome.fasta

    """
}


process Gubbins {

    publishDir "${params.out}/n_${genome}_num_${repeat_range}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_num_${repeat_range}_$filename" }

     input:
        path wholegenome
        each genome
        each repeat_range
        
    
    output:
        path "gubbins.final_tree.tre" , emit: gubbinstree
    
    """   
     run_gubbins  -r GTRGAMMA -p gubbins   ${wholegenome} 
      
    """
}


process get_raxml_tree {

    publishDir "${params.out}/n_${genome}_num_${repeat_range}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_num_${repeat_range}_$filename" }

    input:
        path wholegenome
        each genome
        each repeat_range
    
    output:
        path 'RAxML_bestTree.tree', emit: myRaxML
    
    
    """
     raxmlHPC -m GTRGAMMA   -p 12345 -s ${wholegenome} -N 10 -n tree 
    """
}



process CFML {

    publishDir "${params.out}/n_${genome}_num_${repeat_range}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_num_${repeat_range}_$filename" }

     input:
        path wholegenome
        path myRaxML 
        each genome
        each repeat_range
        
    
    output:
        path "CFML.labelled_tree.newick" , emit: CFMLtree
        path "CFML.importation_status.txt" , emit: CFML_recom
    
    """ 
    
     ClonalFrameML ${myRaxML} ${wholegenome} CFML >  CFML.result.txt  
      
    """
}


process phyloHMM {

     publishDir "${params.out}/n_${genome}_num_${repeat_range}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_num_${repeat_range}_$filename" }

     input:
        path wholegenome
        path myRaxML 
        path recomlog
        path CFML_recom
        each genome
        each nu_hmm
        each mix_prob
        each repeat_range

     output:
        path 'RecomPartial.xml' , emit: partial_XML 
        path 'originalSeq.xml' , emit: original_XML 
        path 'rmse.rmse' , emit : rmse 
        path 'PhyloHMM_Recombination.jpeg' , emit: phyloHMMFig

     
     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/phyloHMM.py -t ${myRaxML} -a ${wholegenome} -nu ${nu_hmm} -x ${params.xml_file} -l ${recomlog} -c ${CFML_recom} -p${mix_prob}
       
        
     """
}


process Beast_alignment {

     publishDir "${params.out}/n_${genome}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_nu_${nu_hmm}_$filename" }

     input:
         path original_XML
         each genome
         each nu_hmm
         
         
     output:    

         path 'wholegenome.trees' , emit:  beastSeqTree
     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/beast  ${original_XML}        
     """
}


process treeannotator_alignment {

     publishDir "${params.out}/n_${genome}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_nu_${nu_hmm}_$filename" }

     input:
     
         path beastSeqTree
         each genome
         each nu_hmm        
         
     output:     
         path 'beastSeqTree.nexus' , emit: SeqTree
     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/treeannotator -b 10  ${beastSeqTree}  beastSeqTree.nexus      
     """
}


process convertor_SeqTree {

     publishDir "${params.out}/n_${genome}" , mode: 'copy' , saveAs: { filename -> "n_${genome}_nu_${nu_hmm}_$filename" }
     
     input:
        path SeqTree
        each genome
        each nu_hmm
    
     output:
         path 'beastSeqTree.newick' , emit : beastTree
         
     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/NexusToNewick.py -t ${SeqTree} -o 'beastSeqTree.newick'
        
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
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/beast  ${partial_XML}        
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
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/treeannotator -b 10  ${beastPartialTree}  beastOurTree.nexus      
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
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/NexusToNewick.py -t ${beastOurTree} -o 'beastHMMTree.newick'
        
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
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/mergeFiles.py ${beastHMMTree}  ${myRaxML}  ${beastTree}  ${CFMLtree}  ${gubbinstree} > allOtherTrees.newick
       
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
       java -jar /home/nehleh/Documents/0_Research/Software/TreeCmp_v2.0-b76/bin/treeCmp.jar  -r ${clonaltree}  -i ${allOtherTrees} -d qt pd rf ms um rfw gdu -o TreeCmpResult.result -W
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
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/rmse_summary.py -f ${rmse}
       
     """
}




workflow {

    BaciSim(genome,repeat_range)
    
    seq_gen(BaciSim.out.BaciSimtrees,frequencies,rates,BaciSim.out.n_genome,BaciSim.out.n_repeat)
    
    Gubbins(seq_gen.out.wholegenome,BaciSim.out.n_genome,BaciSim.out.n_repeat)

    get_raxml_tree(seq_gen.out.wholegenome,BaciSim.out.n_genome,BaciSim.out.n_repeat)
    
    CFML(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,BaciSim.out.n_genome,BaciSim.out.n_repeat)
    
    phyloHMM(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,BaciSim.out.recomlog,CFML.out.CFML_recom,BaciSim.out.n_genome,nu_hmm,mix_prob,BaciSim.out.n_repeat)
    
    Beast_alignment(phyloHMM.out.original_XML,BaciSim.out.n_genome,nu_hmm)
    
    treeannotator_alignment(Beast_alignment.out.beastSeqTree,BaciSim.out.n_genome,nu_hmm)

    convertor_SeqTree(treeannotator_alignment.out.SeqTree,BaciSim.out.n_genome,nu_hmm)   
    
    Beast_partial(phyloHMM.out.partial_XML,BaciSim.out.n_genome,nu_hmm)
    
    treeannotator_partial(Beast_partial.out.beastPartialTree,BaciSim.out.n_genome,nu_hmm)
    
    convertor_ourTree(treeannotator_partial.out.beastOurTree,BaciSim.out.n_genome,nu_hmm)
    
    mergeTreeFiles(convertor_ourTree.out.beastHMMTree,get_raxml_tree.out.myRaxML,convertor_SeqTree.out.beastTree,Gubbins.out.gubbinstree,CFML.out.CFMLtree,BaciSim.out.n_genome,nu_hmm )

    TreeCmp(BaciSim.out.clonaltree,mergeTreeFiles.out.allOtherTrees,BaciSim.out.n_genome,nu_hmm)

    RMSE_summary(BaciSim.out.n_genome,phyloHMM.out.rmse,BaciSim.out.n_repeat)
       
}
