import os
import files_management as fm

def call_hhmake(obj):
    obj.hmm_file = obj.folder+obj.name+".hmm"
    os.system("hhmake -i "+self.train_msa+" -o "+self.hmm_file)

def call_hhfilter(input_file, output_file, hhid):
    print("calling hhfilter "+str(hhid)+" on "+input_file)
    os.system("hhfilter -i "+input_file+" -o "+output_file+" -id "+str(hhid))


def call_trimal(input_file, output_file, trimal_gt, colnumbering_file):
    print("calling trimal "+str(trimal_gt)+" on "+input_file)
    os.system("trimal -in "+input_file+" -out "+output_file+" -gt "+str(trimal_gt)+" -colnumbering > "+colnumbering_file)


def call_reformat(input_file, output_file):
    print("will reformat "+input_file+" using "+fm.get_script_path()+"./reformat.pl to "+output_file)
    call = fm.get_script_path()+"./reformat.pl a3m fas "+input_file+" "+output_file+" -r"
    print(call)
    os.system(call)
