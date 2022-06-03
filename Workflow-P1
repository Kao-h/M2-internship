import glob
from operator import pos
import pandas as pd
import os
import ast
import requests
import gzip
import shutil
from gemmi import cif
import csv
import json

#store Folders
folderPath = "/home/kaouther/RNA_Proj/coconetFinal/"
outputFolder = folderPath + 'JSON_files/'


#extract list of PDB ids from CoConet
pdbIds = pd.read_csv("/home/kaouther/RNA_Proj/coconetFinal/RFam_ID_chain.csv")
pdbIdsList = list(pdbIds['PDB'])
chainList = (list(pdbIds['Chain']))

#download JSON files from DNATCO
def download_JSON(ids,outputFolder):
    """
       This function downloads the corresponding JSON, from a given list of PDB ids, from DNATCO
       and places the output in a given path)
   """
   for pdbid in ids:
      try:
      # download files from url
         url = 'https://dnatco.datmos.org/v4.1/RCSB/json/' + pdbid + '_v41C35A23.json.gz'
         r = requests.get(url, allow_redirects=True)

         # Save the content with name.
         open(outputFolder + pdbid, 'wb').write(r.content)
         with gzip.open(outputFolder + pdbid, 'rb') as f_in:
            with open(outputFolder + pdbid + '.json', 'wb') as f_out:
               shutil.copyfileobj(f_in, f_out)
         os.remove(outputFolder + pdbid)

      except : #if the pdbid is not found by DNATCO
         print(' pass')
         pass

#call function
download_JSON(pdbIdsList, outputFolder)

#create df with data from csv
def createDataset(id_list, outputfolder):
   """
   :param id_list: list of PDB ids
   :param outputfolder: path of output files
   :return: a dataframe

   This function downloads as a CSV file data from DNATCO and puts them in a dataframe format
   extracting the diseried columns
   """

    #columns of bigdf
    columnsName = ["PDB_ID", "chain", "Res1Name", "Res1Pos", "Res2Name", "Res2Pos", "NtCs", "CANA", "nearest_NtC",
                   "pos-pos", "RFams"]

    for pdbid in id_list :

        try:
            # create a df per structure
            rawDataFrame = pd.DataFrame(columns=columnsName)
            # url from DNATCO for output files
            url = 'https://dnatco.datmos.org/v4.1/RCSB/csv/' + pdbid + '_v41C35A23.csv'

            # extract output from dnatco
            rawDataFrame = pd.read_csv(url)
            # new df with columns as in list and index same as raw df
            extracted_df = pd.DataFrame(columns=columnsName, index=rawDataFrame.index)
            # extract data
            extracted_df['PDB_ID'] = rawDataFrame.iloc[0:, 0][0][0:4]
            extracted_df['NtCs'] = rawDataFrame.iloc[0:, 13]
            extracted_df['CANA'] = rawDataFrame.iloc[0:, 14]
            extracted_df["nearest_NtC"] = rawDataFrame.iloc[0:, 18]
            for i in extracted_df.index:
                # extract data from first column
                extracted_df['chain'][i] = rawDataFrame.iloc[0:, 0][i].split('_')[1]
                extracted_df["Res1Name"][i] = rawDataFrame.iloc[0:, 0][i].split('_')[2]
                extracted_df["Res1Pos"][i] = rawDataFrame.iloc[0:, 0][i].split('_')[3]
                extracted_df["Res2Name"][i] = rawDataFrame.iloc[0:, 0][i].split('_')[4]
                extracted_df["Res2Pos"][i] = rawDataFrame.iloc[0:, 0][i].split('_')[5]
                extracted_df["pos-pos"] = extracted_df['Res1Pos'] + '-' + extracted_df['Res2Pos']
            # export as csv new df with only extracted infos
            extracted_df.to_csv(outputfolder + 'df_' + pdbid, index=False, )
        except:  # if the pdbid is not found by DNATCO
            print(' pass')
            pass
    return extracted_df


completeDataset = createDataset(pdbIdsList, folderPath+'/outputdfs/')

# list of merged files returned
def mergeFiles(path, outputfolder):
   """
   :param path: list of PDB ids
   :param outputfolder: path of output files
   :return: a dataframe of merged files

   This function merges files in a given path and concatenate them in a new exported file
   stored in the outputfolder and returns a dataframe
   """
   files =  list(glob.glob(path))
   # joining files with concat and read_csv
   mergedDF = pd.concat(map(pd.read_csv, files), ignore_index=True)  # (85664, 12)
   mergedDF.to_csv(outputfolder, index=False)
   return mergedDF

mergedDF = mergeFiles(folderPath+'/outputdfs/',folderPath+"/NtCsComplete.csv" )

# dict for PDB ids and their corresponding RFAMilies
CoCoNet_rfam_dict = dict(zip(list(pdbIdsList['Rfam']), list(pdbIdsList['PDB'])))
#reverse keys and values in dico
rfam_dict = {v: k for k, v in CoCoNet_rfam_dict.items()}

#add rfams to dataframe
def append_Rfam_in_df(df_initial, dico, existColumn, addColumn):
   """
   :param df_initial : inpurt dataframe
   :param dico: dictionary containing RFAM & PBDs
   :param existColumn : PDB id column
   :param addColumn: column to be added
   :return: a dataframe with the new column

   This function is gonna add a new column to an existing dataframe using a dictionary
   the dictionary keys and the dataframe need to have a common value.
   """

   list_new_column = []
   for row in df_initial[existColumn]:

      temp = None
      for key, values in dico.items():
         if row == values:
            temp = key
      if temp is not None:
         list_new_column.append(temp)
      else:
         list_new_column.append(pd.NA)
   df_initial[addColumn] = list_new_column

   return df_initial


#dowload mmCIF files
def downloadmmCIF(ids,outputFolder)
   """
      :param ids : list of PDB ids
      :param outputFolder: path for output folder
   
   this function downloads mmCIF files for a given list of PDB ids from RCSB
   merges all the files and returns a list of the mmCIF files
   """
   n=0
   for pdbid in ids:
      try:
      # download files from url
         url = 'https://files.rcsb.org/download/' + pdbid + '.cif'
         r = requests.get(url, allow_redirects=True)
         # Save the content with name.
         open(outputFolder + pdbid, 'wb').write(r.content)
      except : #if the pdbid is not found by DNATCO
         print(' pass')
         pass
   completmmCIFfiles = list(glob.glob(outputFolder + '*'))
   return completmmCIFfiles

outputFolder = "/home/kaouther/RNA_Proj/coconetFinal/mmCIF/"
mmCIFfilesList = downloadmmCIF(pdbIdsList, outputFolder)

def mmCIFParser(mmCIFfilesList):
   """:param: mmCIFfilesList list of mmCIF files

   This function parses mmCIF files to extract the list of entities and their corresponding
   sequence for each PDB id
   """
   for mmcifFile in mmCIFfilesList:

      pdbid = mmcifFile.split('/')[6]
      doc = cif.read(mmcifFile)
      block = doc.sole_block()
      listeEntities =list(block.find_values('_entity_poly.entity_id'))
      listeSeq = list(block.find_values('_entity_poly.pdbx_seq_one_letter_code'))
      cleanliste = [None] * len(listeSeq)
      for i in range(len(listeSeq)):
         rem1 = listeSeq[i].replace('\n','')
         rem2 = rem1.replace(';','')
         cleanliste[i] = rem2

      #dict key entity value seq
      dict = {}
      for i in range(len((listeEntities))):
         dict[listeEntities[i]] = cleanliste[i]

      a_file = open("/home/kaouther/RNA_Proj/coconetFinal/mmCIFdict/"+pdbid+"_dict.csv", "w")

      writer = csv.writer(a_file)
      for key, value in dict.items():
         writer.writerow([key, value])
      a_file.close()


#list of dowloaded JSON files

jsonFilesList = list(glob.glob(folderPath + 'JSON_files/'))

def jsonPositions(stepsList, pdbID):
   i = 0
   for step in stepsList:
      part1 = "jq "
      part2 = "'.steps"
      part3 = '."' + step + '".label_seq_id_1,.steps."' + step + '"' + ".label_seq_id_2'" + " " + fileJson
      # print(part1+part2+part3)
      commandInput = part1 + part2 + part3
      tupless = os.popen(commandInput).read().split('\n')
      del tupless[-1]
      # entity
      part4 = '."' + step + '".label_entity_id_1' + "' " + fileJson
      commandInput = part1 + part2 + part4
      entity = os.popen(commandInput).read().split('\n')[0]

      df.loc[i] = [step, step.split('_')[3] + '_' + step.split('_')[5], tupless, str(entity)]
      i += 1
   # print(df.head())
   return df

def JSONfileDF(jsonFilesList):
   """
   :param: jsonFilesList list of dowloaded files
    this function generates dataframe with the steps from DNATCO with the desired columns
   """
   for filename in jsonFilesList:
      # test
      fileJson = filename
      fileName = filename.split('/')[5]
      pdbID = fileName.split('.')[0]
      # Opening JSON file steps
      # commandInput = "jq '.steps' " +fileName
      # data = os.system(commandInput)
      f = open(fileJson)
      # returns JSON object as
      # a dictionary
      data = json.load(f)
      filesSteps = list(data['C5O3_dict'].keys())
      df = pd.DataFrame(columns=['step', 'PDBpos-pos', 'seqPos-seqPos', 'Entity'])
      stepsDataFrame = jsonPositions(filesSteps, fileJson)
      outputfolder = "/home/kaouther/RNA_Proj/coconetFinal/seqPos/"
      stepsDataFrame.to_csv(outputfolder + 'df_' + pdbID.lower(), index=False, )
      print(pdbID)

# correlate json df and mmCIF dict
jsonDFs = list(glob.glob("/home/kaouther/RNA_Proj/coconetFinal/cleanSeq/*"))

def mapJSONfilemmCIF(jsonDFs):
   """
   :param: jsonDFs
   :return: a list of the mapped files list

   This function maps the JSON numbering of the sequence position with their real position in the sequence
   present in the mmCIF file using a common key between the two files
   """
   for file in jsonDFs:

       testfileJson = file
       pdb = testfileJson.split("/")[-1].split(".")[0][0:4]
       testmmCIffile = "/home/kaouther/RNA_Proj/coconetFinal/mmCIFdict/" + pdb + "_dict.csv"
       df = pd.read_csv(testfileJson)
       dico = pd.read_csv(testmmCIffile, index_col=0, header=None, squeeze=True).to_dict()
       EntityColumnJson = list(df.iloc[:, -1])
       EntityList1 = []
       SequenceCol = []

       for element in EntityColumnJson:
           EntityList1.append(int(float(element.replace('"', ''))))

       for element in EntityList1:
           SequenceCol.append(dico[element])

       NtCsdf = pd.read_csv("/home/kaouther/RNA_Proj/coconetFinal/NtCs/df_" + pdb)
       seqPosList = []
       seqList = []
       df['Sequence'] = SequenceCol

       for i in range(len(NtCsdf)):
           step = str(NtCsdf.iloc[i]['PDB_ID']) + '_' + str(NtCsdf.iloc[i]['chain']) + '_' + str(NtCsdf.iloc[i]['Res1Name']) + '_' + str(
               NtCsdf.iloc[i]['Res1Pos']) + '_' + str(NtCsdf.iloc[i]['Res2Name']) + '_' + str(NtCsdf.iloc[i]['Res2Pos'])
           print(step)
           try:
               seqPos = df.loc[df['step'] == step]["seqPos-seqPos"].iloc[0]
               seqPosasList = ast.literal_eval(seqPos)
               realSeqPos = seqPosasList[0] + '-' + seqPosasList[1]
               seqPosList.append(realSeqPos)
               seq = ""
               seq = df.loc[df['step'] == step]["Sequence"].iloc[0]
               if "'" in seq:
                   seq = seq.replace("'", "")
               seqList.append(seq)

           except:
               print(step, 'not found -out ')
       NtCsdf["SeqPos-SeqPos"] = seqPosList
       NtCsdf["Sequence"] = seqList
       NtCsdf["RFams"] = ''

       final_df = append_Rfam_in_df(NtCsdf, CoCoNet_rfam_dict,'PDB_ID', 'RFams')
       final_df.to_csv("/home/kaouther/RNA_Proj/coconetFinal/completedDFs/final_"+pdb, index=False)
       # list of merged files returned
       mappedFilesList = list(glob.glob("/home/kaouther/RNA_Proj/coconetFinal/completedDFs/final_*"))*
       return mappedFilesList

mappedFilesList = mapJSONfilemmCIF(jsonDFs)

#extract only needed column and leave only the corresponding chain in the MSAs
def lightDF(mappedFilesList):
    for file in mappedFilesList:
        df = pd.read_csv(file)
        pdb = file[-4:]
        #drop unneeded columns
        df.drop(['CANA', 'nearest_NtC', 'pos-pos', 'Sequence'], axis=1,inplace=True)
        msa_chain = rfam_dict[pdb]
        df.drop(df[df.chain != msa_chain].index, inplace=True)
        df.to_csv("/home/kaouther/RNA_Proj/coconetFinal/lightDF/light_"+pdb+'.csv', index=False)

lightDF(mappedFilesList)
