import os
import basic_modules.files_management as fm


def call_hhmake(obj):
    obj.hmm_file = obj.folder/(obj.name+".hmm")
    os.system("hhmake -i "+str(self.train_msa)+" -o "+str(self.hmm_file))

def call_hhfilter(input_file, output_file, hhid):
    print("calling hhfilter "+str(hhid)+" on "+str(input_file))
    os.system("hhfilter -i "+str(input_file)+" -o "+str(output_file)+" -id "+str(hhid))


# mieux gérer trimal -> cf Q12404
def call_trimal(input_file, output_file, trimal_gt, cons, colnumbering_file):
    print("calling trimal gt "+str(trimal_gt)+" cons "+str(cons)+" on "+str(input_file))
    os.system("trimal -in "+str(input_file)+" -out "+str(output_file)+" -gt "+str(trimal_gt)+" -cons "+str(cons)+" -colnumbering > "+str(colnumbering_file))


def call_reformat(input_file, output_file):
    print("will reformat "+str(input_file)+" thanks to Soeding's reformat.pl to "+str(output_file))
    call = "reformat.pl a3m fas "+str(input_file)+" "+str(output_file)+" -r"
    print(call)
    os.system(call)
