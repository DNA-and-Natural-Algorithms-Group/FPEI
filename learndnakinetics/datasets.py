import myenums
import learndnakinetics

"""Cisse, Ibrahim I., Hajin Kim, and Taekjip Ha. "A rule of seven in Watson-Crick base-pairing of mismatched sequences." Nature structural & molecular biology 19.6 (2012): 623."""
def read_Cisse2012(done_queue, dataset_list, n , directories, theta_simulation  , load , save  )  :

	for association, j  in  zip([False  ],   ["sup_Fig4g_" ] ):
		bimolecular_reaction = association

		if association == True :
			reaction_type = myenums.ReactionType.HELIXASSOCIATION.value
		else :
			reaction_type = myenums.ReactionType.HELIXDISSOCIATION.value
		dataset_type = myenums.DatasetType.MISMATCH.value

		my_name = '/helix4_cisse2012/'+ j  + str(int(association))
		dataset_path, document , row = learndnakinetics.initconf(my_name , directories )
		temperature_type =myenums.DatasetSpecifications.TEMPCELCIUS.value
		rate_type = myenums.DatasetSpecifications.RATECONSTANT.value
		dataset_name= 	myenums.DatasetName.CISSE2012.value
		learndnakinetics.check_directories (directories)
		docID= dataset_name+ str(bimolecular_reaction)
		dataset= learndnakinetics.Reaction (  counter_celllist =  	[2,3,4,5,6] , num_simulations = 200,  temperature_type =temperature_type  ,   rate_type =  rate_type, dataset_name =dataset_name,   docID=  docID,document =  document, theta_simulation = theta_simulation , reaction_type =  reaction_type, dataset_type= dataset_type , bimolecular_reaction =  bimolecular_reaction , row  = row, dataset_path= dataset_path ,  load = load, save = save   )
		(dataset_list, done_queue, n )  =learndnakinetics.objective_function_auxilary(   done_queue, dataset_list, n ,dataset)
	return (dataset_list, done_queue, n )

"""Hata, Hiroaki, Tetsuro Kitajima, and Akira Suyama. "Influence of thermodynamically unfavorable secondary structures on DNA hybridization kinetics." Nucleic Acids Research (2017)."""
def read_Hata2017 (done_queue,dataset_list, n ,   directories, theta_simulation  ,  load , save ) :

	for association in [True   ]:
		bimolecular_reaction = association
		if association == True :
			reaction_type = myenums.ReactionType.HELIXASSOCIATION.value
			dataset_type = myenums.DatasetType.NODANGLE.value
		else :
			reaction_type = myenums.ReactionType.HELIXDISSOCIATION.value
			dataset_type = myenums.DatasetType.NODANGLE.value
		my_name = '/helix2_hata2017/Table' +str(int(association ))
		dataset_path, document , row = learndnakinetics.initconf(my_name , directories  )
		temperature_type = myenums.DatasetSpecifications.TEMPCELCIUS.value
		rate_type = myenums.DatasetSpecifications.RATECONSTANT10POW5.value
		dataset_name= 	myenums.DatasetName.HATA.value
		learndnakinetics.check_directories (directories)
		docID= dataset_name+ str(bimolecular_reaction)
		dataset= learndnakinetics.Reaction ( counter_celllist = [1, 8, 18 , 23  ], num_simulations = 20, temperature_type =temperature_type  ,   rate_type =  rate_type, dataset_name =dataset_name,   docID=  docID,document =  document, theta_simulation = theta_simulation , reaction_type =  reaction_type, dataset_type= dataset_type , bimolecular_reaction =  bimolecular_reaction , row  = row, dataset_path= dataset_path ,  load = load, save = save )
		(dataset_list, done_queue, n )  =learndnakinetics.objective_function_auxilary(  done_queue, dataset_list, n ,dataset)
	return (dataset_list, done_queue, n )


"""Bonnet, Gragoire, Oleg Krichevsky, and Albert Libchaber. "Kinetics of conformational fluctuations in DNA hairpin-loops." Proceedings of the National Academy of Sciences 95.15 (1998): 8602-8606."""
def read_Bonnet98(done_queue, dataset_list, n , directories, theta_simulation  , load , save )   :

	bimolecular_reaction = False
	for j in [4] :
		for hairpinclosing in [True  , False ]:

			my_name = '/hairpin_bonnet98/Fig'+str(j)  + '_' + str(int(hairpinclosing))
			dataset_path, document , row = learndnakinetics.initconf(my_name , directories)
			learndnakinetics.check_directories (directories)
			dataset_type = myenums.DatasetType.NODANGLE.value
			if  hairpinclosing == True :
				reaction_type = myenums.ReactionType.HAIRPINCLOSING.value
				counter_celllist =  [31,32,33,34,35]
			else :
				counter_celllist = [32,33,34,35,36]
				reaction_type = myenums.ReactionType.HAIRPINOPENING.value
			dataset_name = myenums.DatasetName.BONNET.value
			docID = dataset_name+ str(j) + str(hairpinclosing)
			temperature_type =myenums.DatasetSpecifications.TEMPKELVININV.value
			rate_type =  myenums.DatasetSpecifications.RATECONSTANT.value
			dataset= learndnakinetics.Reaction ( counter_celllist = counter_celllist,  num_simulations =1 ,  temperature_type =temperature_type  ,   rate_type =  rate_type, dataset_name =dataset_name,   docID=  docID,document =  document, theta_simulation = theta_simulation , reaction_type =  reaction_type, dataset_type= dataset_type , bimolecular_reaction =  bimolecular_reaction , row  = row, dataset_path= dataset_path  , load  = load, save = save)
			(dataset_list, done_queue, n )  =learndnakinetics.objective_function_auxilary( done_queue, dataset_list, n ,dataset)
	return (dataset_list, done_queue, n )
