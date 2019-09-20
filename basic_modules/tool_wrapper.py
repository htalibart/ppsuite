import os
import basic_modules.files_management as fm


def call_hhmake(obj):
    obj.hmm_file = obj.folder/(obj.name+".hmm")
    os.system("hhmake -i "+str(self.train_msa)+" -o "+str(self.hmm_file))

def call_hhfilter(input_file, output_file, hhid):
    print("calling hhfilter "+str(hhid)+" on "+str(input_file))
    os.system("hhfilter -i "+str(input_file)+" -o "+str(output_file)+" -id "+str(hhid))


def call_hhblits(input_file, output_file, database, maxfilt=100000, realign_max=100000, B=100000, Z=100000, n=3, e=0.001, retry_hhblits_with_memory_limit_if_fail=False,**kwargs):
    """ calls HH-blits with arguments recommended for CCMpred : https://github.com/soedinglab/CCMpred/wiki/FAQ """
    print("calling hhblits on "+str(input_file)+" using "+str(database)+", output will be available at "+str(output_file))
    hhblits_call = "hhblits -maxfilt "+str(maxfilt)+" -realign_max "+str(realign_max)+" -d "+str(database)+" -all -B "+str(B)+" -Z "+str(Z)+" -n "+str(n)+" -e "+str(e)+" -i "+str(input_file)+" -oa3m "+str(output_file)
    print(hhblits_call)
    os.system(hhblits_call)
    if not output_file.exists() and retry_hhblits_with_memory_limit_if_fail:
        print("HH-blits failed for some reason, trying with a memory limit")
        memory_friendly_call = hhblits_call+" -cpu 1 -maxmem 1"
        os.system(memory_friendly_call)
        if not output_file.exists():
            raise Exception("HH-blits failed. Protein is probably too long ?")


def call_trimal(input_file, output_file, trimal_gt, cons, colnumbering_file):
    print("calling trimal gt "+str(trimal_gt)+" cons "+str(cons)+" on "+str(input_file))
    os.system("trimal -in "+str(input_file)+" -out "+str(output_file)+" -gt "+str(trimal_gt)+" -cons "+str(cons)+" -colnumbering > "+str(colnumbering_file))


def call_reformat(input_file, output_file):
    print("will reformat "+str(input_file)+" thanks to Soeding's reformat.pl to "+str(output_file))
    call = "reformat.pl a3m fas "+str(input_file)+" "+str(output_file)+" -r"
    print(call)
    os.system(call)
