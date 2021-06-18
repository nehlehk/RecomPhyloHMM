nextflow.enable.dsl = 2


params.genomeSize = '5000'
params.recom_len = '700'
params.recom_rate = '0.02'
params.tMRCA = '0.01'
params.nu_sim = '0.2'

params.BaciSimtrees = '/home/nehleh/Desktop/sisters/new_BaciSim/BaciSimTrees.tree'
params.recomlog = '/home/nehleh/Desktop/sisters/new_BaciSim/BaciSim_Log.txt'

params.xml_file = '/home/nehleh/Documents/GTR_template.xml'
params.out = '/home/nehleh/work/results/'


genome = Channel.value(10)
frequencies = Channel.value(' 0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
nu_hmm = Channel.of(0.03)
mix_prob = Channel.of(0.9)
repeat_range = Channel.value(1..1)







process BaciSim {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
     maxForks 1
     
    
     input:
         val genome 
         each repeat_range


     output:
         tuple val(repeat_range), path("BaciSimTrees.tree") , emit: BaciSimtrees
         tuple val(repeat_range), path('clonaltree.tree'), emit: clonaltree
         tuple val(repeat_range), path('BaciSim_Log.txt') , emit: recomlog
         tuple val(repeat_range), path('BaciSim_Recombination.jpeg'), emit: SimFig

         
          
     
     """
       python3.8  /home/nehleh/PhyloCode/RecomPhyloHMM/bin/BaciSim.py -n ${genome} -g ${params.genomeSize} -l ${params.recom_len} -r ${params.recom_rate}  -nu ${params.nu_sim}  -t ${params.tMRCA}
        
     """
}


process seq_gen {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
    maxForks 1   

    input:
        each repeat_range
        val f 
        val r 
          
    output:
        path "wholegenome_${repeat_range}.fasta" , emit: wholegenome
        val repeat_range , emit: range

    
    """ 
     numTrees=\$(wc -l < ${params.BaciSimtrees} | awk '{ print \$1 }')
     seq-gen  -mGTR  -l${params.genomeSize} -r$r -f$f -s0.2 -of ${params.BaciSimtrees} -p \$numTrees > wholegenome_${repeat_range}.fasta

    """
}


process Gubbins {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }

     input:
        path wholegenome
        val repeat_range
        
    
    output:
        path "gubbins.final_tree.tre" , emit: gubbinstree
        
    
    """   
     run_gubbins  -r GTRGAMMA -p gubbins   ${wholegenome} 
      
    """
}


process get_raxml_tree {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
    maxForks 1

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

    input:
        path wholegenome
        path myRaxML
        val repeat_range

     output:
        path 'originalSeq.xml' , emit: original_XML 


     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/MakeXMLSeq.py -t ${myRaxML} -a ${wholegenome}  -x ${params.xml_file}        
        
     """
}


process Beast_alignment {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }


     input:
         path original_XML
         val repeat_range
         
         
     output:    

         path 'wholegenome.trees' , emit:  beastSeqTree
     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/beast  ${original_XML}        
     """
}


process treeannotator_alignment {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }  

     
     input:
     
         path beastSeqTree 
         val repeat_range     
         
     output:     
         path 'beastSeqTree.nexus' , emit: SeqTree
     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/treeannotator -b 10  ${beastSeqTree}  beastSeqTree.nexus      
     """
}


process convertor_SeqTree {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }  

     
     input:
        path SeqTree
        val repeat_range
    
     output:
         path 'beastSeqTree.newick' , emit : beastTree
         
     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/NexusToNewick.py -t ${SeqTree} -o 'beastSeqTree.newick'
        
     """
}



process CFML {

   publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }

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

     input:
        path wholegenome
        path myRaxML 
        path CFML_recom
        path CFMLtree
        val repeat_range


     output:
        path 'CFML_Recombination.jpeg' , emit: CFMLFig
        path 'rmse_CFML.csv' , emit : rmse_CFML 


     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/CMFL_result.py -t ${myRaxML} -a ${wholegenome} -ct ${CFMLtree} -c ${CFML_recom}  -l ${params.recomlog} 
       
        
     """
}


process phyloHMM_two {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }
     maxForks 1

     input:
        path wholegenome
        path myRaxML 
        tuple val(repeat_range), path('recomlog') 
        path CFML_recom
        val repeat_range
        each nu_hmm
        each mix_prob


     output:
        path 'RecomPartial_two.xml' , emit: partial_XML_two
        path 'rmse_phylohmm_two.csv' , emit : rmse_phylohmm_two 
        path 'PhyloHMM_Recombination_two.jpeg' , emit: phyloHMMFig_two
        path 'Recom_phyloHMM_Log_two.txt' , emit: phyloHMMLog_two
        val nu_hmm , emit: my_nu
        val mix_prob , emit: prob
        

     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/phyloHMM_2states.py -t ${myRaxML} -a ${wholegenome}  -x ${params.xml_file} -l ${recomlog}  -c ${CFML_recom} -nu ${nu_hmm} -p ${mix_prob}
        
     """
}




process phyloHMM_four {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }
     maxForks 1

     input:
        path wholegenome
        path myRaxML 
        tuple val(repeat_range), path('recomlog') 
        path CFML_recom
        val repeat_range
        each nu_hmm
        each mix_prob


     output:
        path 'RecomPartial_four.xml' , emit: partial_XML_four 
        path 'rmse_phylohmm_four.csv' , emit : rmse_phylohmm_four
        path 'PhyloHMM_Recombination_four.jpeg' , emit: phyloHMMFig_four
        path 'Recom_phyloHMM_Log_four.txt' , emit: phyloHMMLog_four
        val nu_hmm , emit: my_nu
        val mix_prob , emit: prob
        

     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/phyloHMM.py -t ${myRaxML} -a ${wholegenome}  -x ${params.xml_file} -l ${recomlog}  -c ${CFML_recom} -nu ${nu_hmm} -p ${mix_prob}
        
     """
}



// process RMSE_summary {

//      publishDir "/home/nehleh/work/results/Summary_Results", mode: "copy"
//      maxForks 1
     

//      input: 
//         path collectedRMSE_HMM
//         path collectedRMSE_CFML

  
//      output:
//          path 'RMSE_comparison.jpeg' , emit: rmse_plot
    
//      """
//        python3.8  /home/nehleh/PhyloCode/RecomPhyloHMM/bin/rmse_summary.py -p rmse_phylohmm.csv -c rmse_CFML.csv
       
//      """
// }



process RMSE_summary_states {

     publishDir "/home/nehleh/work/results/Summary_Results", mode: "copy"
     maxForks 1
     

     input: 
        path collectedRMSE_HMM_two
        path collectedRMSE_HMM_four
        path collectedRMSE_CFML

  
     output:
         path 'RMSE_comparison_states.jpeg' , emit: rmse_plot
    
     """
       python3.8  /home/nehleh/PhyloCode/RecomPhyloHMM/bin/rmse_plot_states.py -t rmse_phylohmm_two.csv  -f rmse_phylohmm_four.csv -c rmse_CFML.csv
       
     """
}


process Beast_partial_four {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }


     input:
        path partial_XML_four
        val repeat_range
        each nu_hmm
        each mix_prob
         
         
     output:    
         path 'wholegenome.trees' , emit: beastPartialTree_four
     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/beast  ${partial_XML_four}        
     """
}




process treeannotator_partial_four {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }


     input:   
        path beastPartialTree_four
        val repeat_range
        each nu_hmm
        each mix_prob
         
         
     output:   
     
         path 'beastOurTree_four.nexus' , emit: beastOurTree_four

     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/treeannotator -b 10  ${beastPartialTree_four}  beastOurTree_four.nexus      
     """
}



process convertor_ourTree_four {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }

     
     input: 
        path beastOurTree_four
        val repeat_range
        each nu_hmm
        each mix_prob

    
     output:
         path 'beastHMMTree_four.newick' , emit:  beastHMMTree_four
         
     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/NexusToNewick.py -t ${beastOurTree_four} -o 'beastHMMTree_four.newick'
        
     """
}


process Beast_partial_two {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }


     input:
        path partial_XML_two
        val repeat_range
        each nu_hmm
        each mix_prob
         
         
     output:    
         path 'wholegenome.trees' , emit: beastPartialTree_two
     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/beast  ${partial_XML_two}        
     """
}




process treeannotator_partial_two {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }


     input:   
        path beastPartialTree_two
        val repeat_range
        each nu_hmm
        each mix_prob
         
         
     output:   
     
         path 'beastOurTree_two.nexus' , emit: beastOurTree_two

     
     """
       /home/nehleh/Documents/0_Research/Software/BEAST_with_JRE.v2.6.3.Linux/beast/bin/treeannotator -b 10  ${beastPartialTree_two}  beastOurTree_two.nexus      
     """
}



process convertor_ourTree_two {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }

     
     input: 
        path beastOurTree_two
        val repeat_range
        each nu_hmm
        each mix_prob

    
     output:
         path 'beastHMMTree_two.newick' , emit:  beastHMMTree_two
         
     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/NexusToNewick.py -t ${beastOurTree_two} -o 'beastHMMTree_two.newick'
        
     """
}


process mergeTreeFiles {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }
     maxForks 1

    input:
         path beastHMMTree_four
         path beastHMMTree_two
         path myRaxML 
         path beastTree
         path gubbinstree
         path CFMLtree
         val repeat_range
         each nu_hmm
         each mix_prob    
         

    output:
         path 'allOtherTrees.newick' , emit: allOtherTrees
     
     """
       python3.8 /home/nehleh/PhyloCode/RecomPhyloHMM/bin/mergeFiles.py ${beastHMMTree_four}  ${beastHMMTree_two}  ${myRaxML}  ${beastTree}  ${CFMLtree}  ${gubbinstree} > allOtherTrees.newick
       
     """


}


process TreeCmp {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }
     maxForks 1

     input:
     
         tuple val(repeat_range), path('clonaltree')
         path allOtherTrees
         val repeat_range
         each nu_hmm
         each mix_prob
         
         
     output:
         path 'TreeCmpResult.result' , emit: Comparison     
    
     """
       java -jar /home/nehleh/Documents/0_Research/Software/TreeCmp_v2.0-b76/bin/treeCmp.jar  -r ${clonaltree}  -i ${allOtherTrees} -d qt pd rf ms um rfw gdu -o TreeCmpResult.result -W
     """
}



process TreeCmp_summary {

     publishDir "/home/nehleh/work/results/Summary_Results", mode: "copy"
     maxForks 1
     

     input: 
        path Comparison

  
     output:
         path 'treecmp_Quartet.jpeg' ,       emit: Quartet_plot
         path 'treecmp_PathDiffernce.jpeg' , emit: PathDiffernce_plot
         path 'treecmp_RF.jpeg' ,            emit: RF_plot
         path 'treecmp_MatchingSplit.jpeg' , emit: MatchingSplit_plot
         path 'treecmp_UMAST.jpeg' ,         emit: UMAST_plot
         path 'treecmp_RFWeighted.jpeg' ,    emit: RFWeighted_plot
         path 'treecmp_GeoUnrooted.jpeg' ,   emit: GeoUnrooted_plot
    
     """
       python3.8  /home/nehleh/PhyloCode/RecomPhyloHMM/bin/cmpTree_plot.py -c all_cmpTrees.result
       
     """
}








workflow {

    // BaciSim(genome,repeat_range)
    
    seq_gen(repeat_range,frequencies,rates)
    
    // Gubbins(seq_gen.out.wholegenome,seq_gen.out.range)

    get_raxml_tree(seq_gen.out.wholegenome,seq_gen.out.range)
    
    // make_xml_seq(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,seq_gen.out.range)
    
    // Beast_alignment(make_xml_seq.out.original_XML,seq_gen.out.range)
    
    // treeannotator_alignment(Beast_alignment.out.beastSeqTree,seq_gen.out.range)

    // convertor_SeqTree(treeannotator_alignment.out.SeqTree,seq_gen.out.range)  
   
    CFML(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,seq_gen.out.range)

    CFML_result(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,CFML.out.CFML_recom,CFML.out.CFMLtree,repeat_range)
   
    // phyloHMM_two(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,BaciSim.out.recomlog,CFML.out.CFML_recom,seq_gen.out.range,nu_hmm,mix_prob)
    
    // phyloHMM_four(seq_gen.out.wholegenome,get_raxml_tree.out.myRaxML,BaciSim.out.recomlog,CFML.out.CFML_recom,seq_gen.out.range,nu_hmm,mix_prob)
     
    // collectedRMSE_HMM_two = phyloHMM_two.out.rmse_phylohmm_two.collectFile(name:"rmse_phylohmm_two.csv",storeDir:"/home/nehleh/work/results/Summary_Results", keepHeader:false , sort: false) 
    
    // collectedRMSE_HMM_four = phyloHMM_four.out.rmse_phylohmm_four.collectFile(name:"rmse_phylohmm_four.csv",storeDir:"/home/nehleh/work/results/Summary_Results", keepHeader:false , sort: false)
    
    // collectedRMSE_CFML = CFML_result.out.rmse_CFML.collectFile(name:"rmse_CFML.csv",storeDir:"/home/nehleh/work/results/Summary_Results", keepHeader:false , sort: false) 
    
    // RMSE_summary_states(collectedRMSE_HMM_two,collectedRMSE_HMM_four,collectedRMSE_CFML)
    
    // // RMSE_summary(collectedRMSE_HMM,collectedRMSE_CFML)
    
    
      
    // Beast_partial_four(phyloHMM_four.out.partial_XML_four,seq_gen.out.range,nu_hmm,mix_prob)
    
    // treeannotator_partial_four(Beast_partial_four.out.beastPartialTree_four,seq_gen.out.range,nu_hmm,mix_prob)
    
    // convertor_ourTree_four(treeannotator_partial_four.out.beastOurTree_four,seq_gen.out.range,nu_hmm,mix_prob)
    
    
    // Beast_partial_two(phyloHMM_two.out.partial_XML_two,seq_gen.out.range,nu_hmm,mix_prob)
    
    // treeannotator_partial_two(Beast_partial_two.out.beastPartialTree_two,seq_gen.out.range,nu_hmm,mix_prob)
    
    // convertor_ourTree_two(treeannotator_partial_two.out.beastOurTree_two,seq_gen.out.range,nu_hmm,mix_prob)
    
    
    
    
    // mergeTreeFiles(convertor_ourTree_four.out.beastHMMTree_four,convertor_ourTree_two.out.beastHMMTree_two,get_raxml_tree.out.myRaxML,convertor_SeqTree.out.beastTree,Gubbins.out.gubbinstree,CFML.out.CFMLtree,seq_gen.out.range,nu_hmm,mix_prob)

    // TreeCmp(BaciSim.out.clonaltree,mergeTreeFiles.out.allOtherTrees,seq_gen.out.range,nu_hmm,mix_prob)
    
    // collectedCMP_tree = TreeCmp.out.Comparison.collectFile(name:"all_cmpTrees.result",storeDir:"/home/nehleh/work/results/Summary_Results", keepHeader:false , sort: false) 
    
    // TreeCmp_summary(collectedCMP_tree)
     

          
}