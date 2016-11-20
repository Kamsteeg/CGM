#!/usr/bin/python
#Filename: CGM_Compare_GFF_V101.py
#Part of the compare gene model tool
#Maintainer Jeroen Kamsteeg

# python lib imports
import itertools
import time

class compare_GML:
	def __init__(self, gff3_main_dict, fasta_main_dict, ch_list, MLVL):	
		################ temp vars
		
		# imports
		self.gff3_main_dict = gff3_main_dict
		self.fasta_main_dict = fasta_main_dict
		self.ch_list = ch_list		
		
		# local vars
		self.gm_match_list = []	
		
		# start class functions		
		# -------------------------------------------		
		
		# run the script
		self.compareGML(MLVL)
		
		# -------------------------------------------

	# compare the gff3 gene model items 
	def compareGML(self, MLVL):	
		# go tru the gene model dicts
		for ref, pre in itertools.combinations(self.gff3_main_dict, 2):			
			# loop true the referance gene model dict
			for ref_GM in range(1, len(self.gff3_main_dict[ref])):	
				# range is CH	
				if ref_GM in self.ch_list[ref]:				
					loc = self.ch_list[ref].index(ref_GM)				
					try:
						s_range = [self.ch_list[ref][loc], self.ch_list[ref][loc+1]]	
					except:
						pass		
				# the size of the referance gene model
				ref_start = int(self.gff3_main_dict[ref][ref_GM][0][3])				
				ref_end = int(self.gff3_main_dict[ref][ref_GM][0][4])
				ref_size = ref_end - ref_start				
				# version info
				ref_version = self.gff3_main_dict[ref][ref_GM][0][0]			
				# get the number of exons
				ref_exon = self.gff3_main_dict[ref][ref_GM][1].count("exon") 
				# get the number of CDS 
				ref_cds = self.gff3_main_dict[ref][ref_GM][1].count("CDS")		
				# get the fasta id
				ref_gm_id = self.gff3_main_dict[ref][ref_GM][0][7]	
				
				#get fasta seq 				
				try:	
					ref_gm_seq = self.fasta_main_dict[ref][ref_gm_id].replace('*','')
				except:
					ref_gm_seq = ref
					
				# var
				hit = 0
			
				# loop true the pre ch
				for pre_GM in range(s_range[0], s_range[1]):					
					try:					
						# size of the prediction gene model
						pre_start = int(self.gff3_main_dict[pre][pre_GM][0][3])
						pre_end = int(self.gff3_main_dict[pre][pre_GM][0][4])
						pre_size = pre_end - pre_start						
						# version info
						pre_version = self.gff3_main_dict[pre][pre_GM][0][0]
						# number of items in the gene model
						pre_len = len(self.gff3_main_dict[pre][pre_GM][1])					
						# get the number of exons
						pre_exon = self.gff3_main_dict[pre][pre_GM][1].count("exon")
						# get the number of CDS 
						pre_cds = self.gff3_main_dict[pre][pre_GM][1].count("CDS")						
						# get the fasta id 
						pre_gm_id = self.gff3_main_dict[pre][pre_GM][0][7]
						
						#get fasta seq 
						try:	
							pre_gm_seq = self.fasta_main_dict[pre][pre_gm_id].replace('*','')							
						except:
							pre_gm_seq = pre_gm_id
							
						# match lvl 1
						if int(MLVL) == 1:
							if pre_size == ref_size:	
								# Does the FASTA SEQ match
								if len(set([ref_gm_seq, pre_gm_seq])) == 1:
									# Add to the Hit dict
									item1 = str("gml_"+str(ref)+"_"+str(ref_GM))
									item2 = str("gml_"+str(pre)+"_"+str(pre_GM))		
									match_list = [item1, item2]										
									self.gm_match_list.extend([match_list])																					
									# control if there is a hit
									hit = 1
									# break
									break							 								
						
						# match lvl 2
						if int(MLVL) == 2:
							findFASTA = ref_gm_seq.find(pre_gm_seq)
							if findFASTA == -1:
								pass
							else:
								# Add to the Hit dict
								item1 = str("gml_"+str(ref)+"_"+str(ref_GM))
								item2 = str("gml_"+str(pre)+"_"+str(pre_GM))		
								match_list = [item1, item2]										
								self.gm_match_list.extend([match_list])													
								# control if there is a hit
								hit = 1
								# break
								break							
						# match lvl 3
						if int(MLVL) == 3:
							findFASTA = ref_gm_seq.find(pre_gm_seq)
							if findFASTA == -1:
								pass
							else:
								# Add to the Hit dict
								item1 = str("gml_"+str(ref)+"_"+str(ref_GM))
								item2 = str("gml_"+str(pre)+"_"+str(pre_GM))		
								match_list = [item1, item2]										
								self.gm_match_list.extend([match_list])													
								# control if there is a hit
								hit = 1
								# break							
						# match lvl 4
						if int(MLVL) == 4:
							if pre_size == ref_size:	
								print pre_size, ref_size
								# break
								break 
						
													
					except:
						pass
						self.gm_match_list.extend(["error"])	
						
				# add no match
				if hit == 0:
					item1 = str("gml_"+str(ref)+"_"+str(ref_GM))
					match_list = [item1]
					self.gm_match_list.extend([match_list])	
					
				# control if there is a hit
				hit = 0	
		
	# makes a file for the results
	def result_file_make(self):
		out_file = open(self.NOUTF, "w")
		return out_file
		
	def writeResults(self):
		pass
		
		
		
											
					
	
					
								
								
						
							
																
				

