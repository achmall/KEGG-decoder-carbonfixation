#!/usr/bin/python
def wlp(ko_match): 
	out_data = {'red. acetyl-CoA pathway' : 0, 'acetyl-CoA synthase': 0, 'CO dehydrogenase': 0, 'methyl branch acetogens': 0, 'formate dehydrogenase':0, 'formate-THF ligase': 0, 'methylene-THF dehydrogenase/methenyl-THF cyclohydrolase': 0, 'methylene-THF reductase': 0, 'methyl-THF : corrinoid/iron sulfur protein methyltransferase': 0, 'methyl branch methanogens': 0, 'formylmethanofuran dehydrogenase':0, 'formyl-MFR:tetrahydromethanopterin formyltransferase':0, 'methenyl-tetrahydromethanopterin cyclohydrolase':0, 'methylene-tetrahydromethanopterin dehydrogenase / hydrogenase':0, 'methylene-tetrahydromethanopterin reductase':0, 'tetrahydromethanopterin S-methyltransferase':0}

	acetogen = 0
	methanogen = 0
	
##methyl branch acetogens
#formate dehydrogenase. K05299/K15022: NADP-dependent, K00122/K00123/K00124/K00126/K00127/22515: NAD-dependent (did only take alpha and beta subunits), K22338-41: NAD/Fd bifurcating, K00125/K22516: F420-dependent, K22015: fdh(hydrogenase), K08348: Quinone dependent fdh. There is a lot mor orthologies though so things might be overlooked here.   
	if  (('K05299' in ko_match and 'K15022' in ko_match) 
		or ('K00123' in ko_match and 'K00124' in ko_match) 
		or ('K22338' in ko_match and 'K22339' in ko_match and 'K22340' in ko_match and 'K22341' in ko_match) 
		or ('K00125' in ko_match and 'K22516' in ko_match)
		or ('K22015' in ko_match )
		or ('K08348' in ko_match )):
		acetogen += 1
		out_data['formate dehydrogenase'] = 1
#formate--tetrahydrofolate ligase
	if  ('K01938' in ko_match):
		acetogen += 1
		out_data['formate-THF ligase'] = 1
#methylenetetrahydrofolate dehydrogenase (NADP+) / methenyltetrahydrofolate cyclohydrolase
	if  ('K01491' in ko_match):
		acetogen += 1
		out_data['methylene-THF dehydrogenase/methenyl-THF cyclohydrolase'] = 1
#methylenetetrahydrofolate reductase (NADPH)
	if  ('K00297' in ko_match):
		acetogen += 1
		out_data['methylene-THF reductase'] = 1
#5-methyltetrahydrofolate corrinoid/iron sulfur protein methyltransferase	 							
	if  ('K15023' in ko_match):
		acetogen += 1
		out_data['methyl-THF : corrinoid/iron sulfur protein methyltransferase'] = 1
		
##methyl branch methanogens	
#formylmethanofuran dehydrogenase
	if  ('K00200' in ko_match and 'K00201' in ko_match and 'K00202' in ko_match and 'K00203' in ko_match and 'K11261' in ko_match and 'K00202' in ko_match):
		methanogen += 1
		out_data['formylmethanofuran dehydrogenase'] = 1
#ftr; formylmethanofuran--tetrahydromethanopterin N-formyltransferase
	if  ('K00672' in ko_match):
		methanogen += 1
		out_data['formyl-MFR:tetrahydromethanopterin formyltransferase'] = 1
#mch; methenyltetrahydromethanopterin cyclohydrolase
	if  ('K01499' in ko_match):
		methanogen += 1
		out_data['methenyl-tetrahydromethanopterin cyclohydrolase'] = 1
#hmd; 5,10-methenyltetrahydromethanopterin hydrogenase OR methylene-tetrahydromethanopterin dehydrogenase
	if  ('K13942' in ko_match or 'K00319' in ko_match):
		methanogen += 1
		out_data['methylene-tetrahydromethanopterin dehydrogenase / hydrogenase'] = 1
#mer; 5,10-methylenetetrahydromethanopterin reductase
	if  ('K00320' in ko_match):
		methanogen += 1
		out_data['methylene-tetrahydromethanopterin reductase'] = 1
#mtrABCDEFGH; tetrahydromethanopterin S-methyltransferase
	if  ('K00577' in ko_match and 'K00578' in ko_match and 'K00579' in ko_match and 'K00580' in ko_match and 'K00581' in ko_match and 'K00582' in ko_match and 'K00583' in ko_match and 'K00584' in ko_match):
		methanogen += 1
		out_data['tetrahydromethanopterin S-methyltransferase'] = 1
					
	
	out_data['methyl branch acetogens'] = (acetogen / 5)* 0.49 + 0.5 #This is to show them in green / show complete subpathways in blue	
	out_data['methyl branch methanogens'] = (methanogen / 6)* 0.49 + 0.5
		
			
##Carbon fixing branch, CODH/ACS
#acetyl-CoA decarbonylase/synthase complex subunit alpha OR CO-methylating acetyl-CoA synthase (14138 is acetogens, 00192 is methanogens) 
	if  ('K00192' in ko_match or 'K14138' in ko_match):
		out_data['acetyl-CoA synthase'] = 1
		acetogen += 1
		methanogen += 1
##catalytic subunits only of CO dehydrogenase
#anaerobic carbon-monoxide dehydrogenase OR aerobic carbon-monoxide dehydrogenase large subunit
	if  ('K00198' in ko_match or 'K03520' in ko_match):
		
		out_data['CO dehydrogenase'] = 1	
		acetogen += 1
		methanogen += 1
		
	acetogencompleteness = (acetogen / 7)* 0.49	
	methanogencompleteness = (methanogen / 8)* 0.49	
	
	if acetogencompleteness > methanogencompleteness:
		out_data['red. acetyl-CoA pathway'] = acetogencompleteness
	else:	
		out_data['red. acetyl-CoA pathway'] = methanogencompleteness
		
	if  (out_data['acetyl-CoA synthase']  == 1 and out_data['CO dehydrogenase'] == 1):
		out_data['red. acetyl-CoA pathway'] += 0.5
			
	return out_data
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
def cbbc(ko_match):
	out_data = {'CBB Cycle': 0,'rubisco':0, 'phosphoglycerate kinase':0, 'glyceraldehyde 3-phosphate dehydrogenase':0,'phosphoribulokinase':0,'fbp/sbp aldolase':0,'fbpase':0,'fbpap':0,'transketolase':0,'sbpase':0,'transaldolase':0,'ribose 5-phosphate isomerase':0,'ribulose-phosphate 3-epimerase':0,'triosephosphate isomerase':0}
	
	cbbc = 0
	
#phosphoglycerate kinase
	if   ('K00927' in ko_match):
		out_data['phosphoglycerate kinase']  = 1
		cbbc += 1
#glyceraldehyde 3-phosphate dehydrogenase 
	if  ('K00134' in ko_match or 'K05298' in ko_match or 'K00150' in ko_match or 'K10705' in ko_match) :
		out_data['glyceraldehyde 3-phosphate dehydrogenase']  = 1
		cbbc += 1
#Ribulose regeneration
#fbp/sbp aldolase + fbpase OR fbpap
	if   ('K01622' in ko_match or 'K01623' in ko_match or 'K01624' in ko_match or 'K11645' in ko_match or 'K16305' in ko_match or 'K16306' in ko_match): 
		out_data['fbp/sbp aldolase']  = 1
	if 	('K11532' in ko_match or 'K03841' in ko_match or 'K02446' in ko_match or 'K01086' in ko_match or 'K04041' in ko_match):
		 out_data['fbpase']  = 1
	if 	('K01622' in ko_match):
		out_data['fbpap']  = 1
	if 	((out_data['fbp/sbp aldolase']  == 1 and out_data['fbpase']  == 1) or (out_data['fbpap']  == 1)):
		cbbc += 1
#transketolase 
	if  ('K00615' in ko_match):
		out_data['transketolase']  = 1
		cbbc += 1
#SBPase OR transaldolase
	if  ('K01086' in ko_match or 'K01100' in ko_match or 'K11532' in ko_match or 'K22315' in ko_match):
		out_data['sbpase']  = 1
	if  ('K00616' in ko_match or 'K13810' in ko_match):
		out_data['transaldolase']  = 1
	if  (out_data['transaldolase']  == 1 or out_data['sbpase']  == 1):
		cbbc += 1
# ribose 5-phosphate isomerase OR ribulose-phosphate 3-epimerase
	if  ('K01807' in ko_match or 'K01808' in ko_match):
		out_data['ribose 5-phosphate isomerase']  = 1
	if 	('K01783' in ko_match):
		out_data['ribulose-phosphate 3-epimerase']  = 1
	if  (out_data['ribose 5-phosphate isomerase']  == 1 or out_data['ribulose-phosphate 3-epimerase'] == 1):
		cbbc += 1
#triosephosphate isomerase
	if  ('K01803' in ko_match):
		out_data['triosephosphate isomerase']  = 1
		cbbc += 1
#RuBisCO - Only large subunit
	if  ('K01601' in ko_match):
		out_data['rubisco'] = 1
		cbbc += 1
#phosphoribulokinase            
	if  ('K00855' in ko_match):
		out_data['phosphoribulokinase']  = 1
		cbbc += 1
	
	out_data['CBB Cycle'] = (float(cbbc) / float (9))* float(0.49) 
		
	if  (out_data['rubisco']  == 1 and out_data['phosphoribulokinase'] == 1):
		out_data['CBB Cycle'] += 0.5
		
		
	return out_data

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------		
def acetylCoAsuccinylCoApws(ko_match): 
	out_data = {'roTCA Cycle' : 0, 'rTCA Cycle' : 0, '3OHP Cycle':0, '4OHB/3OHP Cycle':0, 'DC/4OHB Cycle': 0, 
	
	'citrate branch': 0, 'citrate synthase':0, 'ATP-citrate lyase' : 0,'citryl-CoA synthetase' :0,'citryl-CoA lyase' :0, 'aconitate hydratase' :0, 'isocitrate dehydrogenase':0, '2-oxoglutarate carboxylase':0, '2-oxoglutarate synthase':0, 
	
	'dicarboxylate branch': 0,'pyruvate synthase': 0, 'pyruvate->oxaloacetate':0,'malate dehydrogenase':0,'fumarate hydratase':0,'fumarate reductase':0,'succinyl-CoA synthesis':0,
	
	'3OHP branch': 0, 'bifunctional malonyl-CoA/malonate semialdehyde reductase':0,'malonyl-CoA reductase':0,'malonate semialdehyde reductase':0,'Propionyl-CoA synthase':0,'3OHP-CoA synthetase':0, '3OHP-CoA dehydratase':0, 'acryloyl-CoA reductase':0, 'Acetyl-CoA carboxylase/Propionyl-CoA carboxylase':0, 'methylmalonyl-CoA epimerase':0, 'methylmalonyl-CoA mutase':0,
	
	'Glyoxylate branch': 0,  'succinyl-CoA:(S)-malate CoA-transferase':0, 'succinate dehydrogenase':0, 'fumarate hydratase':0, 'malyl-CoA/(S)-citramalyl-CoA lyase':0, 'mesaconyl-C1-CoA hydratase':0, 'mesaconyl-CoA C1-C4 CoA transferase':0, 'mesaconyl-C4-CoA hydratase':0,
	
	'4OHB branch': 0, 'succinyl-CoA reductase':0, 'succinate semialdehyde reductase':0, '4-OHB-CoA synthetase':0, '4-OHB-CoA dehydratase':0, 'crotonyl-CoA hydratase':0, '3-OHB-CoA dehydrogenase':0, 'acetoacetyl-CoA ß-ketothiolase':0
	
	}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
##citrate branch

	citrate = 0
		
#citrate synthase 
	if ('K01647' in ko_match):
		out_data['citrate synthase'] = 1
#ATP-citrate lyase	
	if ('K15230' in ko_match and 'K15231' in ko_match):
		out_data['ATP-citrate lyase'] = 1
#citryl-CoA synthetase		
	if ('K15232' in ko_match and 'K15233' in ko_match):
		out_data['citryl-CoA synthetase'] = 1
#citryl-CoA lyase
	if	('K15234' in ko_match or 'K01644' in ko_match):	
		out_data['citryl-CoA lyase'] = 1
## citrate cleavate combined	
	if (out_data['ATP-citrate lyase'] == 1 or out_data['citrate synthase'] == 1 or (out_data['citryl-CoA synthetase'] == 1 and out_data['citryl-CoA lyase'] == 1)):
		citrate += 1
		
#aconitate hydratase
	if ('K01681' in ko_match or 'K01682' in ko_match):
		out_data['aconitate hydratase'] = 1
		citrate += 1
#isocitrate dehydrogenase (K00031 includes the H. thermophilus oxalosuccinate reductase)
	if ('K00031' in ko_match or 'K00030' in ko_match or 'K17753' in ko_match):
		out_data['isocitrate dehydrogenase'] = 1
		citrate += 1
#2-oxoglutarate carboxylase (hydrogenobacter enzyme), because we cannot differentiate between oxalosuccinate reductase and isocitrate dehydrogenase this is just for information and doesnt contribute to the score.
	if ('K20140' in ko_match):
		out_data['2-oxoglutarate carboxylase'] = 1	
		
#2-oxoglutarate/2-oxoacid ferredoxin oxidoreductase
	if ('K00174' in ko_match and 'K00175' in ko_match):
		out_data['2-oxoglutarate synthase'] = 1
		citrate += 1
#citrate branch
	citratecompleteness = (citrate / 4)* 0.49 
	out_data['citrate branch'] = citratecompleteness
#key enzymes citrate branch	
	if ((out_data['ATP-citrate lyase'] == 1 or out_data['citrate synthase'] == 1 or (out_data['citryl-CoA synthetase'] == 1 and out_data['citryl-CoA lyase'] == 1))
		and
		(out_data['2-oxoglutarate synthase'] == 1)):
		
		out_data['citrate branch'] += 0.5		
	
	
	
	
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
##Dicarboxylate branch

	dicarboxylate = 0	
	 
#pyruvate ferredoxin oxidoreductase (pyruvate synthase) alpha and beta subunits. More subunits exist and there are many different pyruvate synthases characterized, but all seem to have the 2 biggest subunits, that is why I included these. K03737 is a homodimer enzyme from desulfovibrio.
	if (('K00169' in ko_match and 'K00170' in ko_match)
		or
		('K00174' in ko_match and 'K00175' in ko_match)
		or
		('K03737' in ko_match)):
		out_data['pyruvate synthase'] = 1
		dicarboxylate += 1
#Pyruvate carboxylase OR ((PEP synthase (Pyruvate:water dikinase) or pyruvate, orthophosphate dikinase) and (PEP carboxylase or PEP carboxykinase [it is not sure how well PEPCK functions in oxaloacetate synthesis, but could work])). Left out malic enzyme, because it favours decarboxylation a lot. 
	if (('K01958' in ko_match or ('K01959' in ko_match and 'K01960' in ko_match)) #PC
		or
		(('K01007' in ko_match or 'K01006' in ko_match) #Pyruvate:water dikinase, pyruvate:Pi dikinase
		and
		('K01595' in ko_match or ('K01596' in ko_match or 'K20370' in ko_match or 'K01610' in ko_match)))): #PEPC OR (PEPCK (GTP) OR PEPCK (diphosphate) OR PEPCK (ATP))
		out_data['pyruvate->oxaloacetate'] = 1
		dicarboxylate += 1	
#malate dehydrogenase
	if ('K00116' in ko_match or 'K00024' in ko_match or 'K00025' in ko_match or 'K00026' in ko_match or 'K00051' in ko_match): #quinone, NADH, NADH, NADH, NADPH, To Do: No orthology for archaeal malate dehydrogenase, build hmm from https://www.genome.jp/dbget-bin/www_bget?ec:1.1.1.299 --> https://epub.uni-regensburg.de/32702/: Igni_1263, annotated as 00024, although with lower score . See https://www.genome.jp/entry/iho:Igni_1263. 
		out_data['malate dehydrogenase'] = 1
		dicarboxylate += 1
#fumarate hydratase
	if ('K01676' in ko_match or 'K01677' in ko_match or 'K01678' in ko_match or 'K01679' in ko_match):
		out_data['fumarate hydratase'] = 1
		dicarboxylate += 1
#fumarate reductase 
	#if (('K00244' in ko_match and 'K00245' in ko_match and 'K00246' in ko_match and 'K00247' in ko_match) #quinone
	#	or
	#	('K00239' in ko_match and 'K00240' in ko_match and 'K00241' in ko_match and 'K00242' in ko_match) #quinone
	#	or
	#	('K00234' in ko_match and 'K00235' in ko_match and 'K00236' in ko_match and 'K00237' in ko_match) #SDH (no clear distinction between sdh and fdr, so included this as well)
	#	or
	#	('K18556' in ko_match and 'K18557' in ko_match and 'K18558' in ko_match and 'K18559' in ko_match and 'K18560' in ko_match) #Hydrogenobacter NADH-Dependent frd
	#	or 
	#	('K18209' in ko_match and 'K18210' in ko_match)): #Fumarate reductase (thiol)	2.7.2.8	 	
	#		out_data['fumarate reductase'] = 0
	## I rewrote this to make more sense: The flavoprotein subunit must be there in any case, plus either a menaquinone binding subunit or a nad binding subunit or a heterodisulfide reductase homologous subunit. Omit anchor proteins and iron sulfur proteins because they are very small/maybe not very specific.
	### Removed everything but the flavoprotein subunits --> the other subunits are too variable and are not always annotated correctly by kofamscan
	
	if ( 'K00244' in ko_match or 'K00239' in ko_match or 'K00234' in ko_match or  'K18209' in ko_match):   #flavoprotein 
		#and 
		#(('K00241' in ko_match or 'K00247' in ko_match or 'K00236' in ko_match)  #cytochrome (menaquinone-binding) subunit. FrdD. 
		#or
		#('K18559' in ko_match) #nad binding subunit
		#or
		#('K18210' in ko_match))): # heterodisulfide like
		out_data['fumarate reductase'] = 1
		dicarboxylate += 1
					
#succinyl-CoA synthetase OR acetyl-CoA:succinate CoA-transferase
	if (('K01899' in ko_match and 'K01900' in ko_match) or #GTP dependent
		('K01902' in ko_match and 'K01903' in ko_match) or #ATP dependent
		('K18118' in ko_match)): #acetyl-CoA:succinate CoA-transferase
		out_data['succinyl-CoA synthesis'] = 1
		dicarboxylate += 1
		
	dicarboxylatecompleteness = (dicarboxylate / 6)* 0.49 
	out_data['dicarboxylate branch'] = dicarboxylatecompleteness + 0.5			
	
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
## 3OHP Branch		

	hydroxypropionate = 0

# bifunctional malonyl-CoA/malonate semialdehyde reductase OR malonyl-CoA reductase and malonate semialdehyde reductase K14468
	# bifunctional malonyl-CoA/malonate semialdehyde reductase
	if('K14468' in ko_match):
		out_data['bifunctional malonyl-CoA/malonate semialdehyde reductase'] = 1
	# malonyl-CoA reductase
	if('K15017' in ko_match):
		out_data['malonyl-CoA reductase'] = 1
	# malonate semialdehyde reductase
	if ('K15039' in ko_match or 'K18602' in ko_match or 'K16066' in ko_match or 'K09019' in ko_match):
		out_data['malonate semialdehyde reductase'] = 1
	
	if(out_data['bifunctional malonyl-CoA/malonate semialdehyde reductase'] == 1 or (out_data['malonyl-CoA reductase'] == 1 and out_data['malonate semialdehyde reductase'] == 1 )):
		hydroxypropionate += 1											
# Propionyl-CoA synthase OR 3OHP-CoA synthetase AND 3OHP-CoA dehydratase AND acryloyl-CoA reductase
	# Propionyl-CoA synthase
	if('K14469' in ko_match):
		out_data['Propionyl-CoA synthase'] = 1	
	# 3OHP-CoA synthetase 
	if('K15018' in ko_match or 'K18594' in ko_match):
		out_data['3OHP-CoA synthetase'] = 1
	# 3OHP-CoA dehydratase
	if('K15019' in ko_match):
		out_data['3OHP-CoA dehydratase'] = 1
	# acryloyl-CoA reductase
	if('K15020' in ko_match):
		out_data['acryloyl-CoA reductase'] = 1			
	if(out_data['Propionyl-CoA synthase'] == 1 or (out_data['3OHP-CoA synthetase'] == 1 and out_data['3OHP-CoA dehydratase'] == 1 and out_data['acryloyl-CoA reductase'] == 1)):
		hydroxypropionate += 1													
# Acetyl-CoA carboxylase/Propionyl-CoA carboxylase (chloroflexus/metallosphaera/nmar enzymes)	
	if(('K02160' in ko_match and 'K01961' in ko_match and 'K01962' in ko_match and 'K01963' in ko_match) #these are the chloroflexus enzymes
		or
		('K01964' in ko_match and 'K15037' in ko_match and 'K15036' in ko_match) # these are the Metallosphaera enzymes
		or
		('K18604' in ko_match and 'K18603' in ko_match and 'K18605' in ko_match)): # these are the Nmar enzymes
		out_data['Acetyl-CoA carboxylase/Propionyl-CoA carboxylase'] = 1	
		hydroxypropionate += 1
# methylmalonyl-CoA epimerase		
	if('K05606' in ko_match):
		out_data['methylmalonyl-CoA epimerase'] = 1	
		hydroxypropionate += 1
# methylmalonyl-CoA mutase		
	if(('K01847' in ko_match)
		or
		('K01848' in ko_match and 'K01849' in ko_match)):
		out_data['methylmalonyl-CoA mutase'] = 1	
		hydroxypropionate += 1
	
	hydroxypropionatecompleteness = (hydroxypropionate / 5)* 0.49	
	out_data['3OHP branch'] = hydroxypropionatecompleteness
#key enzymes 3OHP branch --> Removed malonyl-CoA reductase from key enzymes because it is hard to identify and not clear which gene it is in thaumarchaeota
	if (out_data['Acetyl-CoA carboxylase/Propionyl-CoA carboxylase'] == 1 
		#and
		#(out_data['bifunctional malonyl-CoA/malonate semialdehyde reductase'] == 1 or (out_data['malonyl-CoA reductase'] == 1 and out_data['malonate semialdehyde reductase'] == 1 ))
		and
		(out_data['Propionyl-CoA synthase'] == 1 or (out_data['3OHP-CoA synthetase'] == 1 and out_data['3OHP-CoA dehydratase'] == 1 and out_data['acryloyl-CoA reductase'] == 1))
		and
		out_data['methylmalonyl-CoA mutase'] == 1):

		out_data['3OHP branch'] += 0.5
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
##Glyoxylate branch

	glyoxylate = 0

#succinyl-CoA:(S)-malate CoA-transferase
	if('K14471' in ko_match and 'K14472' in ko_match):
		out_data['succinyl-CoA:(S)-malate CoA-transferase'] = 1		
		glyoxylate += 1		
#succinate dehydrogenase 
##removed everything but the flavoproteins, same like for fumarate reductase

	if('K00244' in ko_match or 'K00239' in ko_match or 'K00234' in ko_match):

	#if (('K00244' in ko_match and 'K00245' in ko_match and 'K00246' in ko_match and 'K00247' in ko_match)
	#	or
	#	('K00239' in ko_match and 'K00240' in ko_match and 'K00241' in ko_match and 'K00242' in ko_match)
	#	or
	#	('K00234' in ko_match and 'K00235' in ko_match and 'K00236' in ko_match and 'K00237' in ko_match)): 
		out_data['succinate dehydrogenase'] = 1
		glyoxylate += 1
#fumarate hydratase --> the original keggdecoder entry is wrong because the ecoli fumarase genes each code for an independent fumarase. Put all of the KOs in OR
	if ('K01677' in ko_match or 'K01678' in ko_match or 'K01679' in ko_match or 'K01676' in ko_match):
		out_data['fumarate hydratase'] = 1
		glyoxylate += 1
#trifunctional (S)-malyl-CoA /ß-methylmalyl-CoA /(S)- citramalyl-CoA lyase
	if('K08691' in ko_match): # This is the main enzyme from chloroflexus. There are also other KOs from the methylaspartate cycle/itaconate degradation: K18292 K19281 K19282. Decided not to include them.												
		out_data['malyl-CoA/(S)-citramalyl-CoA lyase'] = 1
		glyoxylate += 1
# mesaconyl-C1-CoA hydratase	4.1.2.13 / 3.1.3.14 (=2-methylfumaryl-CoA hydratase)
	if('K14449' in ko_match):
		out_data['mesaconyl-C1-CoA hydratase'] = 1
		glyoxylate += 1
# mesaconyl-CoA C1-C4 CoA transferase	5.3.1.9 (=2-methylfumaryl-CoA isomerase) 
	if('K14470' in ko_match):
		out_data['mesaconyl-CoA C1-C4 CoA transferase'] = 1
		glyoxylate += 1															
# mesaconyl-C4-CoA hydratase	5.1.3.4 (=3-methylfumaryl-CoA hydratase)														
	if('K09709' in ko_match):
		out_data['mesaconyl-C4-CoA hydratase'] = 1
		glyoxylate += 1
		
	glyoxylatecompleteness = (glyoxylate / 7) * 0.49
	out_data['Glyoxylate branch'] = glyoxylatecompleteness
# key enzyme glyoxylate branch	
	if out_data['malyl-CoA/(S)-citramalyl-CoA lyase'] == 1 :
		out_data['Glyoxylate branch'] += 0.5
	
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
##4OHB branch

	hydroxybutyrate = 0

#succinyl-CoA reductase 	
	if('K15017' in ko_match or 'K15038' in ko_match or 'K18601' in ko_match): #msed OR pyrobaculum OR nmar. Nmar_1608 --> K18601 (Walker, 2010; not confirmed yet). No candidate gene in Ignicoccus. Not NADP(H) dependent (description paper)
		hydroxybutyrate += 1
		out_data['succinyl-CoA reductase'] = 1
#succinate semialdehyde reductase To Do: no orthology for nmar, for igni : https://www.kegg.jp/dbget-bin/www_bget?iho:Igni_0132 (confirmed in Jennifer Flechslers PhD thesis, https://epub.uni-regensburg.de/32702/) annotated as K11173
	if('K14465' in ko_match or 'K11173' in ko_match): #(NADPH) 
		hydroxybutyrate += 1
		out_data['succinate semialdehyde reductase'] = 1
#4-OHB-CoA synthetase    
	if('K14466' in ko_match or 'K14467' in ko_match or 'K18861' in ko_match or 'K18593' in ko_match or 'K01895' in ko_match): # K18593 is ADP forming, K01895 is annotated as an acCoA synthetase.To Do: no orthology for Igni, no clear candidates...Kegg annotates this Igni_0475 as K18861, but in Jennys thesis it wasnt active
		hydroxybutyrate += 1
		out_data['4-OHB-CoA synthetase'] = 1
#4-OHB-CoA dehydratase
	if('K14534' in ko_match):
		hydroxybutyrate += 1
		out_data['4-OHB-CoA dehydratase'] = 1
#crotonyl-CoA hydratase
	if('K15016' in ko_match or 'K15019' in ko_match ): #K15016 is bifunctional crotonase/3-hydroxybutyryl-CoA dehydrogenase, K15019 is the bifunctional 3-OHP-CoA dehydratase/crotonyl-CoA hydratase (Liu, 2021, msphere). To Do: No KOs for Msed_0336 and Msed_0384 (unifunctional crotonases from (Liu, 2020, frontiers)
		hydroxybutyrate += 1
		out_data['crotonyl-CoA hydratase'] = 1
#3-OHB-CoA dehydrogenase. There are a lot more orthologies that could be added here .... see R06941, R01975, R03026. 
	if('K15016' in ko_match or 'K00074' in ko_match):
		hydroxybutyrate += 1
		out_data['3-OHB-CoA dehydrogenase'] = 1
#acetoacetyl-CoA ß-ketothiolase	
	if('K00626' in ko_match):
		hydroxybutyrate += 1
		out_data['acetoacetyl-CoA ß-ketothiolase'] = 1
		
	hydroxybutyratecompleteness = (hydroxybutyrate / 7)* 0.49
	out_data['4OHB branch'] = hydroxybutyratecompleteness
#key enzymes 4OHB branch
		
	if (out_data['4-OHB-CoA dehydratase'] == 1):
		out_data['4OHB branch'] += 0.5
	
				
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
# r(o)TCA cycles
	out_data['roTCA Cycle'] = out_data['rTCA Cycle'] = (citratecompleteness + dicarboxylatecompleteness) / 2 # divided by two, because the two branches combined shoud be max 0.49 when they are complete.
	
#key enzymes rTCA cycle
	if (out_data['2-oxoglutarate synthase'] == 1 
		and 
		(out_data['ATP-citrate lyase'] == 1 or (out_data['citryl-CoA synthetase'] == 1 and out_data['citryl-CoA lyase'] == 1))):
		
		out_data['rTCA Cycle'] += 0.5
	
#key enzymes roTCA Cycle
	if (out_data['2-oxoglutarate synthase'] == 1 
		and 
		out_data['citrate synthase'] == 1):
		
		out_data['roTCA Cycle'] += 0.5
	
#3OHBP cycle
	out_data['3OHP Cycle'] = (hydroxypropionatecompleteness + glyoxylatecompleteness) / 2
		
#key enzymes 3OHP cycle	
	if (out_data['Acetyl-CoA carboxylase/Propionyl-CoA carboxylase'] == 1 
		and
		(out_data['bifunctional malonyl-CoA/malonate semialdehyde reductase'] == 1 or (out_data['malonyl-CoA reductase'] == 1 and out_data['malonate semialdehyde reductase'] == 1 ))
		and
		(out_data['Propionyl-CoA synthase'] == 1 or (out_data['3OHP-CoA synthetase'] == 1 and out_data['3OHP-CoA dehydratase'] == 1 and out_data['acryloyl-CoA reductase'] == 1))
		and
		out_data['methylmalonyl-CoA mutase'] == 1
		and
		out_data['malyl-CoA/(S)-citramalyl-CoA lyase'] == 1):
		
		out_data['3OHP Cycle'] += 0.5
			
#3OHP/4OHB cycle
	out_data['4OHB/3OHP Cycle'] = (hydroxypropionatecompleteness + hydroxybutyratecompleteness) / 2
	
#key enzymes 3OHP/4OHB cycle
	if (out_data['Acetyl-CoA carboxylase/Propionyl-CoA carboxylase'] == 1 
		#and
		#(out_data['bifunctional malonyl-CoA/malonate semialdehyde reductase'] == 1 or (out_data['malonyl-CoA reductase'] == 1 and out_data['malonate semialdehyde reductase'] == 1 ))
		and
		(out_data['Propionyl-CoA synthase'] == 1 or (out_data['3OHP-CoA synthetase'] == 1 and out_data['3OHP-CoA dehydratase'] == 1 and out_data['acryloyl-CoA reductase'] == 1))
		and
		out_data['methylmalonyl-CoA mutase'] == 1
		and
		out_data['4-OHB-CoA dehydratase'] == 1):
	
		out_data['4OHB/3OHP Cycle'] += 0.5

#DC/4OHB cycle
	out_data['DC/4OHB Cycle'] = (dicarboxylatecompleteness + hydroxybutyratecompleteness) / 2
	
#DC/4OHB cycle key enzyme		
	if out_data['4-OHB-CoA dehydratase'] == 1:
		
		out_data['DC/4OHB Cycle'] += 0.5	#key enzyme 4OHb branch
	

	return out_data
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
def redglycpw(ko_match):
	out_data = {'formate dehydrogenase': 0,'formate-THF ligase':0,'methylene-THF dehydrogenase/methenyl-THF cyclohydrolase':0, 'Glycine cleavage/synthase system': 0,'Glycine reductase complex':0,'Glycine reductase complex':0,'acetyl-P --> acetyl-CoA':0,'pyruvate synthase':0,'serine hydroxymethyltransferase':0,'Serine deaminase':0,'Red. Glycine Pathway':0}
	
	glycinereductase = 0
	glycineserine = 0
	

##enzymes from supplementary data 3 from Sánchez-Andrea, 2020 --> used the enzymes from ATCC_27774 to search for corresponding KOs

	if (('K05299' in ko_match and 'K15022' in ko_match) 
		or ('K00123' in ko_match and 'K00124' in ko_match) 
		or ('K22338' in ko_match and 'K22339' in ko_match and 'K22340' in ko_match and 'K22341' in ko_match) 
		or ('K00125' in ko_match and 'K22516' in ko_match)
		or ('K22015' in ko_match )
		or ('K08348' in ko_match )):
		
		out_data['formate dehydrogenase'] = 1
		glycinereductase += 1
		glycineserine += 1
#Formate-THF ligase
#formate--tetrahydrofolate ligase. The Desulfovibrio enzyme (Ddes_1296) doesnt have a KO. used the one thats in kegg. 
	if ('K01938' in ko_match):
		
		out_data['formate-THF ligase'] = 1
		glycinereductase += 1
		glycineserine += 1
#methylenetetrahydrofolate dehydrogenase (NADP+) / methenyltetrahydrofolate cyclohydrolase
	if ('K01491' in ko_match):
		
		out_data['methylene-THF dehydrogenase/methenyl-THF cyclohydrolase'] = 1
		glycinereductase += 1
		glycineserine += 1
#Glycine cleavage/synthase system (dihydrolipoamide dehydrogenase K00382, or L-Protein is not mentioned in the paper. Leave out for now. Conversely, this K0800 is not mentioned in kegg pathway(?)
# --> Paper: The GCS was assumed to be present if three out of four genes could be detected
	if (('K00283' in ko_match and 'K00282' in ko_match) #glycine dehydrogenase (decarboxylating), subunit 2/1, glycine cleavage system P-protein
		and ('K02437' in ko_match) # Glycine cleavage system H-protein
		and ('K00605' in ko_match) # Aminomethyltransferase, glycine cleavage system T-protein
		and ('K03800' in ko_match)): # Lipoyltransferase/lipoate-protein ligase
		
		out_data['Glycine cleavage/synthase system'] = 1
		glycinereductase += 1
		glycineserine += 1
		
#Glycine reductase complex (omitted the thioredoxin genes (K03671 and K00384), see explanation in the paper
# --> Paper: for the GR five out of seven genes (under the assumption that the two thioredoxin genes, or any other two genes, could be present at another genomic location), 
	if (('K10670' in ko_match) 					# Glycine reductase complex, component A
		and ('K10671' in ko_match and 'K10672' in ko_match) 	# Glycine reductase complex, component B, subunits alpha and beta / subunit gamma
		and ('K21576' in ko_match and 'K21577' in ko_match)): 	# Glycine reductase complex, component C, subunit alpha / beta
		
		out_data['Glycine reductase complex'] = 1
		glycinereductase += 1
		
#Acetyl-P --> acetyl-CoA
	if (('K00625' in ko_match or 'K13788' in ko_match or 'K15024' in ko_match) #Phosphate acetyl transferase (K13788 is the one from the paper corresponding to Ddes_0297, added the other ones for 2.3.1.8 from kegg) 
		or
		('K00925' in ko_match) and (('K01895' in ko_match)or('K01913' in ko_match))): #Acetate kinase (only one KO) + #Acetyl-CoA synthetase	
		
		out_data['acetyl-P --> acetyl-CoA'] = 1
		glycinereductase += 1
		
#pyruvate ferredoxin oxidoreductase (pyruvate synthase) alpha and beta subunits. More subunits exist and there are many different pyruvate synthases characterized, but all seem to have the 2 biggest subunits, so lets just go for these ones. K03737 is the desulfovibrio enzyme (Ddes_0298). 
	if (('K00169' in ko_match and 'K00170' in ko_match)
		or
		('K00174' in ko_match and 'K00175' in ko_match)
		or
		('K03737' in ko_match)):
			
		out_data['pyruvate synthase'] = 1
		glycinereductase += 1
			
#serine hydroxymethyltransferase (Ddes_1617)
	if ('K00600' in ko_match):
			out_data['serine hydroxymethyltransferase'] = 1
			glycineserine += 1
#Serine deaminase (NO KO! TO DO build custom hmm profile for Ddes_2155/Serine dehydratase-like/Uncharacterised protein family UPF0597). Used KOs for EC 4.3.1.19   instead.
	if (('K01754' in ko_match)
		or
		('K17989' in ko_match )):
			out_data['Serine deaminase'] = 1
			glycineserine += 1
			
# pathway definition (cited from the paper):
# The pathway was considered to be present 
# if any of the four FDHs was present, 
# together with the GCS, the DsvG11_3068 formate-THF ligase, 
# the DsvG11_1518 methenyl-THF cyclohydrolase 
# and either each of the two optional routes: 
# (A) GR and at least one of phosphate acetyltransferase (DsvG11_0941) or acetate kinase (DsvG11_0942), 
# or (B) glycine hydroxymethyltransferase (DsvG11_2276) and SDA (DsvG11_1577). 

	glycinereductasecompleteness = (glycinereductase / 7) * 0.49
	glycineserinecompleteness = (glycineserine / 6) * 0.49

	if (glycinereductasecompleteness > glycineserinecompleteness):
		
		out_data['Red. Glycine Pathway'] = glycinereductasecompleteness
	else:	
		
		out_data['Red. Glycine Pathway'] = glycineserinecompleteness

#key enzymes
			
	if (out_data['Glycine cleavage/synthase system']==1
		and
		(out_data['Glycine reductase complex']==1 or (out_data['serine hydroxymethyltransferase']==1 and out_data['Serine deaminase']==1))):
			
		out_data['Red. Glycine Pathway'] += 0.5	
	
	return out_data


##accessory proteins that were upregulated in the paper, didnt include
#Nitrogen regulatory protein P-II		Ddes_1148|nitrogen regulatory protein P-II
#Glutamate-ammonia ligase			Ddes_1149|glutamine synthetase catalytic region
#Glutamate synthase, small subunit		Ddes_1157|glutamate synthase alpha subunit domain protein
#Glutamate synthase, large subunit		Ddes_1156|ferredoxin-dependent glutamate synthase
#Glutamine amidotransferase type 2 like	Ddes_1155|putative glutamate synthase, large subunit
#glutamate dehydrogenase (NADP+)		Ddes_1152|Glutamate dehydrogenase (NADP(+))
#formate/nitrite transporter			Ddes_0933|formate/nitrite transporter K06212
#Ammonium transporter				Ddes_1147|ammonium transporter K03320
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
def default_viz(genome_df, outfile_name):
	import seaborn as sns
	import matplotlib
	import matplotlib.pyplot as plt
#------------define my custom colour scale--------------------
	cdict = {'red': ((0.0, 1.0, 1.0),
	                 (0.5, 0.8, 1.0),   
	                 (0.99, 0.0, 0.3),
	                 (1.0, 0.3, 0.3)),
	         'green': ((0.0, 1.0, 1.0),
	                  (0.5, 0.2, 1.0),
	                  (0.99, 0.8, 0.3),
	                  (1.0, 0.3, 0.3)),
	         'blue': ((0.0, 1.0, 1.0),
	                  (0.5, 0.2, 1.0),
	                  (0.99, 0.0, 1.0),
	                  (1.0, 1.0, 1.0))}
	my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,1000)
#-------------------------------------------------------------	
	sns.set(font_scale=1.2)
	sns.set_style({"savefig.dpi": 200})
	ax = sns.heatmap(genome_df, cmap=my_cmap, linewidths=2,
		linecolor='k', square=True, xticklabels=True,
		yticklabels=True, cbar_kws={"shrink": 0.1})
	ax.xaxis.tick_top()
	#ax.set_yticklabels(ax.get_yticklabels(), rotation=90)
	plt.xticks(rotation=90)
	plt.yticks(rotation=0)
	# get figure (usually obtained via "fig,ax=plt.subplots()" with matplotlib)
	fig = ax.get_figure()
	# specify dimensions and save
	#xLen = len(genome_df.columns.values.tolist())*20
	#yLen = len(genome_df.index.tolist())*20
	fig.set_size_inches(100, 100)
	fig.savefig(outfile_name, bbox_inches='tight', pad_inches=0.1)

def main():
	import os
	import matplotlib
	matplotlib.use('Agg')
	import argparse
	import pandas as pd
	from scipy.cluster import hierarchy
	from scipy.spatial import distance


	parser = argparse.ArgumentParser(description="Accepts KEGG KOALA\
									text file as input. Produces function\
									list and heat map figure.")
	parser.add_argument('-i', '--input', help="Input KOALA file. See documentation\
						for correct format")
	parser.add_argument('-t', '--tangleopt', help="Number of tree iterations for minimizing tangles in tanglegram", default=1000)
	parser.add_argument('-o', '--output', help="List version of the final heat\
						map figure")
	parser.add_argument('-v', '--vizoption', help="Options: static, interactive, tanglegram")
	parser.add_argument('--newick', help="Required input for tanglegram visualization")
	parser.add_argument("-m", "--myorder", help ="Orders output as specified by	user.", default="None")
	parser.add_argument("-p", "--pathway", help ="Options: all, detail, wjp, cbbc, accoasuccoa, redglyc", default="all")
	args = parser.parse_args()
	arg_dict = vars(args)

	genome_data = {}

	for line in open(str(arg_dict['input']), "r"):
		line = line.rstrip()
		info = line.split()
		if len(info) > 1:
			if info[0].split("_")[0] in genome_data.keys():
				genome_data[info[0].split("_")[0]].append(info[1])
			else:
				genome_data[info[0].split("_")[0]] = [info[1]]
#All Lines in the heatmap have to be entered here, exactly like the out_data in the pathway definitions, and in the right order
	if arg_dict['pathway'] == 'all':
		function_order = ['red. acetyl-CoA pathway','CBB Cycle','roTCA Cycle','rTCA Cycle','3OHP Cycle','4OHB/3OHP Cycle','DC/4OHB Cycle','Red. Glycine Pathway']
	if arg_dict['pathway'] == 'detail':
		function_order = ['red. acetyl-CoA pathway','acetyl-CoA synthase','CO dehydrogenase','methyl branch acetogens','formate dehydrogenase','formate-THF ligase','methylene-THF dehydrogenase/methenyl-THF cyclohydrolase','methylene-THF reductase','methyl-THF : corrinoid/iron sulfur protein methyltransferase','methyl branch methanogens','formylmethanofuran dehydrogenase','formyl-MFR:tetrahydromethanopterin formyltransferase','methenyl-tetrahydromethanopterin cyclohydrolase','methylene-tetrahydromethanopterin dehydrogenase / hydrogenase','methylene-tetrahydromethanopterin reductase','tetrahydromethanopterin S-methyltransferase','CBB Cycle','rubisco','phosphoglycerate kinase','glyceraldehyde 3-phosphate dehydrogenase','phosphoribulokinase','fbp/sbp aldolase','fbpase','fbpap','transketolase','sbpase','transaldolase','ribose 5-phosphate isomerase','ribulose-phosphate 3-epimerase','triosephosphate isomerase','roTCA Cycle','rTCA Cycle','3OHP Cycle','4OHB/3OHP Cycle','DC/4OHB Cycle','citrate branch','citrate synthase','ATP-citrate lyase','citryl-CoA synthetase','citryl-CoA lyase','aconitate hydratase','isocitrate dehydrogenase','2-oxoglutarate carboxylase','2-oxoglutarate synthase','dicarboxylate branch','pyruvate synthase','pyruvate->oxaloacetate','malate dehydrogenase','fumarate hydratase','fumarate reductase','succinyl-CoA synthesis','3OHP branch','bifunctional malonyl-CoA/malonate semialdehyde reductase','malonyl-CoA reductase','malonate semialdehyde reductase','Propionyl-CoA synthase','3OHP-CoA synthetase','3OHP-CoA dehydratase','acryloyl-CoA reductase','Acetyl-CoA carboxylase/Propionyl-CoA carboxylase','methylmalonyl-CoA epimerase','methylmalonyl-CoA mutase','Glyoxylate branch','succinyl-CoA:(S)-malate CoA-transferase','succinate dehydrogenase','malyl-CoA/(S)-citramalyl-CoA lyase','mesaconyl-C1-CoA hydratase','mesaconyl-CoA C1-C4 CoA transferase','mesaconyl-C4-CoA hydratase','4OHB branch','succinyl-CoA reductase','succinate semialdehyde reductase','4-OHB-CoA synthetase','4-OHB-CoA dehydratase','crotonyl-CoA hydratase','3-OHB-CoA dehydrogenase','acetoacetyl-CoA ß-ketothiolase','Red. Glycine Pathway','Glycine cleavage/synthase system','Glycine reductase complex','acetyl-P --> acetyl-CoA','serine hydroxymethyltransferase','Serine deaminase']
	if arg_dict['pathway'] == 'wjp':
		function_order = ['red. acetyl-CoA pathway','acetyl-CoA synthase','CO dehydrogenase','methyl branch acetogens','formate dehydrogenase','formate-THF ligase','methylene-THF dehydrogenase/methenyl-THF cyclohydrolase','methylene-THF reductase','methyl-THF : corrinoid/iron sulfur protein methyltransferase','methyl branch methanogens','formylmethanofuran dehydrogenase','formyl-MFR:tetrahydromethanopterin formyltransferase','methenyl-tetrahydromethanopterin cyclohydrolase','methylene-tetrahydromethanopterin dehydrogenase / hydrogenase','methylene-tetrahydromethanopterin reductase','tetrahydromethanopterin S-methyltransferase']
	if arg_dict['pathway'] == 'cbbc':
		function_order = ['CBB Cycle','rubisco','phosphoglycerate kinase','glyceraldehyde 3-phosphate dehydrogenase','phosphoribulokinase','fbp/sbp aldolase','fbpase','fbpap','transketolase','sbpase','transaldolase','ribose 5-phosphate isomerase','ribulose-phosphate 3-epimerase','triosephosphate isomerase']
	if arg_dict['pathway'] == 'accoasuccoa':
		function_order = ['roTCA Cycle','rTCA Cycle','3OHP Cycle','4OHB/3OHP Cycle','DC/4OHB Cycle','citrate branch','citrate synthase','ATP-citrate lyase','citryl-CoA synthetase','citryl-CoA lyase','aconitate hydratase','isocitrate dehydrogenase','2-oxoglutarate carboxylase','2-oxoglutarate synthase','dicarboxylate branch','pyruvate synthase','pyruvate->oxaloacetate','malate dehydrogenase','fumarate hydratase','fumarate reductase','succinyl-CoA synthesis','3OHP branch','bifunctional malonyl-CoA/malonate semialdehyde reductase','malonyl-CoA reductase','malonate semialdehyde reductase','Propionyl-CoA synthase','3OHP-CoA synthetase','3OHP-CoA dehydratase','acryloyl-CoA reductase','Acetyl-CoA carboxylase/Propionyl-CoA carboxylase','methylmalonyl-CoA epimerase','methylmalonyl-CoA mutase','Glyoxylate branch','succinyl-CoA:(S)-malate CoA-transferase','succinate dehydrogenase','malyl-CoA/(S)-citramalyl-CoA lyase','mesaconyl-C1-CoA hydratase','mesaconyl-CoA C1-C4 CoA transferase','mesaconyl-C4-CoA hydratase','4OHB branch','succinyl-CoA reductase','succinate semialdehyde reductase','4-OHB-CoA synthetase','4-OHB-CoA dehydratase','crotonyl-CoA hydratase','3-OHB-CoA dehydrogenase','acetoacetyl-CoA ß-ketothiolase']
	if arg_dict['pathway'] == 'redglyc':
		function_order = ['Red. Glycine Pathway','Glycine cleavage/synthase system','Glycine reductase complex','acetyl-P --> acetyl-CoA','serine hydroxymethyltransferase','Serine deaminase']


	filehandle = str(arg_dict['output'])
	out_file = open(filehandle, "w")
	out_file.write('Function'+"\t"+str("\t".join(function_order))+"\n")
#All pathway functions have to be entered here
	for k in genome_data:
		pathway_data = {}
		pathway_data.update(wlp(genome_data[k]))
		pathway_data.update(cbbc(genome_data[k]))
		pathway_data.update(acetylCoAsuccinylCoApws(genome_data[k]))
		pathway_data.update(redglycpw(genome_data[k]))
	#    print k, pathway_data

		out_string = str(k)+"\t"
		out_list = [k]
		for i in function_order:
			out_list.append(pathway_data[i])
		out_string = str(out_list).strip('[]')
		tab_string = ""
		for l in out_string:
			if l == "\'":
				continue
			if l == ",":
				tab_string = tab_string + "\t"
			else:
				tab_string = tab_string + l
		out_file.write(tab_string+"\n")
	out_file.close()


	file_in = open(filehandle, "r")
	genome = pd.read_csv(file_in, index_col=0, sep='\t')
	rearrange = False
	if arg_dict["myorder"] != 'None' and os.path.exists(arg_dict["myorder"]):
		rearrange = True
		leaf_order = []
		for line in open(str(arg_dict["myorder"]), "r"):
			line = line.rstrip("\r\n")
			leaf_order.append(line)
		genome = genome.reindex(leaf_order)

	if arg_dict['vizoption'] == 'static':
		from .KEGG_clustering import hClust_euclidean
		#from KEGG_clustering import hClust_euclidean
		if len(genome.index) >= 2 and not rearrange:
			genome = hClust_euclidean(genome)
		default_viz(genome, os.path.splitext(filehandle)[0] + ".svg")
	if arg_dict['vizoption'] == 'interactive':
		from .Plotly_viz import plotly_viz
		#from Plotly_viz import plotly_viz
		plotly_viz(genome, os.path.splitext(filehandle)[0] + ".html")
	if arg_dict['vizoption'] == 'tanglegram':
		from .MakeTanglegram import make_tanglegram
		#from MakeTanglegram import make_tanglegram
		if len(genome.index) >= 3:
			make_tanglegram(genome, str(arg_dict['newick']), os.path.splitext(filehandle)[0] + ".tanglegram.svg", int(arg_dict["tangleopt"]))
		else:
			raise ValueError("Tanglegram mode requires three or more genomes")


if __name__ == "__main__":
	main()
