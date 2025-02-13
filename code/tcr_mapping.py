import editdistance
import os
import pandas as pd
from tqdm import tqdm
import pickle as pkl

def get_barcode_position(bead_bc_file,bead_pos_file):
    def low_cpx_bc(bc,max_homo):
        return Counter(list(bc)).most_common(1)[0][1]>max_homo
    tmp = open(bead_pos_file).readlines()
    pos_x = [float(it) for it in tmp[0].strip().split(",")]
    pos_y =  [float(it) for it in tmp[1].strip().split(",")]
    bc_list = [it.strip().replace(",","") for it in open(bead_bc_file).readlines()]
    max_homo = int(0.8*len(list(bc_list[0])))  # assume all have the same length
    bc_pos_dict = {}
    for bc,coordx,coordy in zip(bc_list,pos_x,pos_y):
        if not low_cpx_bc(bc,max_homo):
            bc_pos_dict[bc] = (coordx,coordy)
    return bc_pos_dict, bc_list

def build_6mer_dist(bc_list):
    start_km = {}
    mid_km = {}
    end_km = {}
    for bc in bc_list:
        start_km.setdefault(bc[:6] , []).append(bc)
        mid_km.setdefault(bc[4:10], []).append(bc)
        end_km.setdefault(bc[-6:] , []).append(bc)
    return start_km,mid_km,end_km

def barcode_matching(bc_pos_dict,spatial_bc_list,max_dist=2):
    bc_matching_dict = {}
    def get_sel_bc(bc):
        res = []
        if bc[:6] in start_km:
            res += start_km[bc[:6]]
        if bc[-6:] in end_km:
            res += end_km[bc[-6:]]
        if bc[4:10] in mid_km:
            res += mid_km[bc[4:10]]
        return set(res)
    exact_match = 0
    fuzzy_match =0
    bc_ref_list = list(bc_pos_dict.keys())
    start_km,mid_km,end_km = build_6mer_dist(bc_ref_list)
    for bc in spatial_bc_list:
        if bc in bc_pos_dict:
            exact_match += 1
            bc_matching_dict[bc] = bc
        else:
            sel_bc = get_sel_bc(bc)
            if len(sel_bc)>0:
                fz = [(it, editdistance.eval(it, bc)) for it in sel_bc]
                fz = [it for it in fz if it[1]<=max_dist]
                fz.sort(key=lambda x:x[1])
                if len(fz)==0:
                    continue
                if len(fz)>1 and fz[0][1]==fz[1][1]:
                    if editdistance.eval(fz[0][0][:-1], bc[-1])>editdistance.eval(fz[1][0][:-1], bc[-1]):  # higher error rate in the last base of the barcode
                        fuzzy_match += 1
                        bc_matching_dict[bc] = fz[1][0]
                    elif editdistance.eval(fz[0][0][:-1], bc[-1])<editdistance.eval(fz[1][0][:-1], bc[-1]):
                        fuzzy_match += 1
                        bc_matching_dict[bc] = fz[0][0]
                else:
                    fuzzy_match += 1
                    bc_matching_dict[bc] = fz[0][0]
    return bc_matching_dict,exact_match,fuzzy_match

class tcr_mapping:
    def __init__(self,trb_folder,puck_folder,sn):
        self.trb_folder = trb_folder
        self.puck_folder = puck_folder
        
        self.abbrev_puck_name = self.puck_folder[11:25]
        self.trb_file_name = f'{self.abbrev_puck_name}_B_result_TRB.tsv'
        self.reads_folder = f'{self.abbrev_puck_name}_B_reads'

        self.bc_matching_folder = f'{self.puck_folder}/barcode_matching'
        self.bc_matching = f'{self.puck_folder[11:]}_barcode_matching.txt'
        self.bc_matching_gz = f'{self.puck_folder[11:]}_barcode_matching.txt.gz'
        self.out_fn = f'{self.abbrev_puck_name}_B'
        self.sample_name =f'{self.abbrev_puck_name}_B'
        self.sn = sn
        self.fastq_file_name = f'{self.abbrev_puck_name}_B_S{self.sn}_R1_001.fastq.gz'
    def find_sample_num(self):
        client = storage.Client()
        BUCKET_NAME = 'fc-secure-2423a594-f5fe-40c6-b998-f28ac17ea012'
        bucket = client.get_bucket(BUCKET_NAME)

        blobs = bucket.list_blobs()

        for blob in blobs:
            if self.trb_folder in blob.name:
                if 'R2_001.fastq' in blob.name:
                    if self.abbrev_puck_name in blob.name:
                        test = blob.name
                        if 'B_S' in test:
                            idx = test.index('B_S')
                            self.sn = test[idx+3]
                        
       
    def load_files(self):
        os.system(f'/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp gs://fc-secure-2423a594-f5fe-40c6-b998-f28ac17ea012/{self.trb_folder}/{self.trb_file_name}  ./')
        os.system(f'/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp gs://fc-secure-2423a594-f5fe-40c6-b998-f28ac17ea012/{self.trb_folder}/{self.fastq_file_name} ./')
        os.system(f'/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil -m -q cp -r -n  gs://fc-secure-2423a594-f5fe-40c6-b998-f28ac17ea012/{self.trb_folder}/{self.reads_folder}/ ./')
        os.system(f'/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil -m cp -r gs://fc-secure-2423a594-f5fe-40c6-b998-f28ac17ea012/{self.bc_matching_folder}/{self.bc_matching} ./')
        os.system(f'/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil -m cp -r gs://fc-secure-2423a594-f5fe-40c6-b998-f28ac17ea012/{self.bc_matching_folder}/{self.bc_matching_gz} ./')
        os.system(f'gunzip -f ./{self.reads_folder}/*.gz')
        os.system(f'gunzip -f ./{self.bc_matching_gz}')
        os.system(f'gunzip -f ./{self.fastq_file_name}') 
    
    def link_reads_to_seq(self):
        self.bc_loc_df = pd.read_csv(self.bc_matching,sep='\t',header=None)
        self.locs = [(i,j) for i, j in zip(self.bc_loc_df[2],self.bc_loc_df[3])]

        self.bc_loc_dict = {i:j for i,j in zip(self.bc_loc_df[0],self.locs)}

        # Identify cell barcodes
        self.bcs = self.bc_loc_df[0]
        self.bcs = {i:1 for i in self.bcs} # Just so it's a dictionary for Luyi's code


        # Build a dictionary linking readID to barcode from sequencing data
        self.fastq_file_r1 = self.fastq_file_name[:-3]

        line_count = 0
        all_read_count = 0
        self.readIDtobarcode = {}

        with open(self.fastq_file_r1,"r") as file1:
            for line in file1:
                if line_count % 4 == 0:  # readID name
                    line = line.strip()
                    read_name = line[1:].split(' ')[0]
                if line_count % 4 == 1:  # Check if fastq sequence
                    all_read_count += 1
                    barcode = line[0:8] + line[26:33]

                    umi = line[33:41]

                    self.readIDtobarcode[read_name] = [barcode,umi]
                line_count += 1
                if line_count % 1000000 == 0:
                    print(line_count)
        file1.close()
        print("Total reads:", all_read_count)
    
    def hamming_correct(self,rerun=True):
        if rerun:
            # Hamming correct the barcode assignments to only look at ones that match the reads
            self.readIDtobarcode_filt = {}
            miss = 0
            read_bcs = []
            for r in tqdm(self.readIDtobarcode):
                bc = self.readIDtobarcode[r][0]
                if bc in self.bcs:
                    self.readIDtobarcode_filt[r] = self.readIDtobarcode[r]
                else: 
                    read_bcs.append(bc)

            read_bcs = set(read_bcs)
            # Try hamming correcting the barcode
            (mapping_dict,_,_) = barcode_matching(self.bcs,read_bcs)

            for r in tqdm(self.readIDtobarcode):
                bc = self.readIDtobarcode[r][0]
                if bc in mapping_dict:
                    self.readIDtobarcode_filt[r] = [mapping_dict[bc] , self.readIDtobarcode[r][1]]
                elif r in self.readIDtobarcode_filt:
                    continue
                else: 
                    miss += 1

            print(miss, 'reads do not have a cell barcode that matches the barcode matching list within 1 hamming distance. ')
            print(len(mapping_dict), 'reads were corrected to a barcode')

            # Save the barcode correction, because it's slow.
            fn = open(f"{self.sample_name}_readIDtobarcode_filt.pkl","wb")
            pkl.dump(self.readIDtobarcode_filt,fn)
            fn.close()
        else:
            with open(f"{self.sample_name}_readIDtobarcode_filt.pkl", 'rb') as f:
                self.readIDtobarcode_filt = pkl.load(f)

    def make_read_name_to_clone(self,rerun=True):
        use_normal = True
        
        
        self.clone_df = pd.read_csv(self.trb_file_name,sep='\t')
        self.cloneIdtoCDR3 = {i:j for i,j in zip(self.clone_df.cloneId,self.clone_df.aaSeqCDR3)}
        self.readIDtocloneID = {}

        cloneId = self.clone_df.cloneId

        if rerun:
            line_count = 0
            for i in tqdm(cloneId):
                do_gunzip = True
                file_name = self.out_fn+'_reads_cln'+str(i)+'.fastq.gz' # R2?
                file_name_mod = self.out_fn+'_reads_cln'+str(i)+'_R2.fastq.gz' # R2?
                
                
                direc = self.out_fn+'_reads'
                fastq_file_r2 = "./"+direc+"/"+ file_name
                fastq_file_r2_mod = "./"+direc+"/"+ file_name_mod
                
                if ((os.path.exists(fastq_file_r2[:-3])) or (os.path.exists(fastq_file_r2))):
                    use_normal = True
                    if os.path.exists(fastq_file_r2[:-3]):
                        do_gunzip = False
                    if do_gunzip:
                        print('file_here')
                        os.system('gunzip -f {fastq_file_r2}')
                        do_gunzip = False
                    
                if ((os.path.exists(fastq_file_r2_mod[:-3])) or (os.path.exists(fastq_file_r2_mod))):
                    use_normal = False
                    if os.path.exists(fastq_file_r2_mod[:-3]):
                        use_normal = False
                        do_gunzip = False
                    if do_gunzip:
                        use_normal = False
                        print('file_here')
                        os.system('gunzip -f {fastq_file_r2_mod}')
                        do_gunzip = False

                if use_normal:
                    if not os.path.exists(fastq_file_r2[:-3]):
                        print(f'{fastq_file_r2[:-3]} file_missing')
                        continue
                else:
                    if not os.path.exists(fastq_file_r2_mod[:-3]):
                        print(f'{fastq_file_r2_mod[:-3]} file_missing')
                        continue
                        
                if use_normal:
                    fastq_file_r2 = fastq_file_r2
                else:
                    fastq_file_r2 = fastq_file_r2_mod
                    
                with open(fastq_file_r2[:-3],"r") as file1:
                    for line in file1:
                        if line_count % 4 == 0:  # readID name
                            line = line.strip()
                            read_name = line[1:].split(' ')[0]
                            self.readIDtocloneID[read_name] = i
                        line_count += 1

            #Save the barcode correction, because it's slow.
            fn = open(f"{self.sample_name}_readIDtoclone.pkl","wb")
            pkl.dump(self.readIDtocloneID,fn)
            fn.close()
            
        else:
            with open(f"{self.sample_name}_readIDtoclone.pkl", 'rb') as f:
                self.readIDtocloneID = pkl.load(f)
                        
    def save_mapped_locs(self):
        self.barcode_to_cloneID = {}

        for read in self.readIDtobarcode_filt:
            barcode = self.readIDtobarcode_filt[read][0]
            if barcode not in self.barcode_to_cloneID:
                self.barcode_to_cloneID[barcode] = []
            if read in self.readIDtocloneID:
                cloneID = self.readIDtocloneID[read]
                self.barcode_to_cloneID[barcode].append(cloneID)
        self.barcode_to_cloneID = {i:list(set(self.barcode_to_cloneID[i])) for i in self.barcode_to_cloneID}

        bcs_with_cdr3 = list(self.barcode_to_cloneID.keys())

        colnames = list(self.clone_df.columns[1:])
        reporter_df = {i:[] for i in colnames}

        x_loc = []
        y_loc = []
    
        barcode = []
        cdr3 = []
        for bc in tqdm(self.barcode_to_cloneID):
            if len(self.barcode_to_cloneID[bc]) != 0:
                for cloneID in self.barcode_to_cloneID[bc]:
                    clone = self.cloneIdtoCDR3[cloneID]

                    barcode.append(bc)
                    x_loc.append(self.bc_loc_dict[bc][0])
                    y_loc.append(self.bc_loc_dict[bc][1])
                    cdr3.append(clone)
                    for c in colnames:
                        reporter_df[c].append(self.clone_df[self.clone_df.cloneId==cloneID][c].values[0])

        final_dict = {'x':x_loc,'y':y_loc,'barcode':barcode,'cdr3':cdr3}

        for c in colnames:
            final_dict[c] = reporter_df[c]

        pd.DataFrame.from_dict(final_dict).to_csv(self.out_fn+'.csv')

        saved_fn = self.out_fn + '.csv'

        os.system(f'/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp ./{saved_fn} gs://fc-secure-2423a594-f5fe-40c6-b998-f28ac17ea012/{self.puck_folder}')
        print('Completed!')
        
        

def pkl_dump(obj,file_name):
    file = open(f'{file_name}.pkl', 'wb')
    # dump information to that file
    pkl.dump(obj, file)
    # close the file
    file.close()  
    
def pkl_load(file_name):
    file = open(f'{file_name}.pkl','rb')
    obj = pkl.load(file)
    file.close()
    return obj