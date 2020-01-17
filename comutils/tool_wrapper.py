import os
import pathlib
import subprocess

from Bio.Align.Applications import MuscleCommandline

import comutils.files_management as fm

def call_hhmake(obj):
    obj.hmm_file = obj.folder/(obj.name+".hmm")
    cmd = "hhmake -i "+str(self.train_msa)+" -o "+str(self.hmm_file)
    subprocess.Popen(cmd, shell=True).wait()

def call_hhfilter(input_file, output_file, hhid):
    print("calling hhfilter "+str(hhid)+" on "+str(input_file))
    fm.check_if_file_ok(input_file)
    cmd = "hhfilter -i "+str(input_file)+" -o "+str(output_file)+" -id "+str(hhid)
    subprocess.Popen(cmd, shell=True).wait()


def call_hhblits(input_file, output_file, database, maxfilt=100000, realign_max=100000, B=100000, Z=100000, n=3, e=0.001, retry_hhblits_with_memory_limit_if_fail=False, hhr_file=None, **kwargs):
    """ calls HH-blits with arguments recommended for CCMpred : https://github.com/soedinglab/CCMpred/wiki/FAQ """
    if hhr_file is None:
        hhr_file = pathlib.Path('.'.join(str(output_file).split('.')[:-1])+".hhr")
    print("calling hhblits on "+str(input_file)+" using "+str(database)+", output will be available at "+str(output_file))
    fm.check_if_file_ok(input_file)
    hhblits_call = "hhblits -maxfilt "+str(maxfilt)+" -realign_max "+str(realign_max)+" -d "+str(database)+" -all -B "+str(B)+" -Z "+str(Z)+" -n "+str(n)+" -e "+str(e)+" -i "+str(input_file)+" -oa3m "+str(output_file)+" -o "+str(hhr_file)
    print(hhblits_call)
    subprocess.Popen(hhblits_call, shell=True).wait()
    if not output_file.exists() and retry_hhblits_with_memory_limit_if_fail:
        print("HH-blits failed for some reason, trying with a memory limit")
        memory_friendly_call = hhblits_call+" -cpu 1 -maxmem 1"
        subprocess.Popen(memory_friendly_call, shell=True).wait()
        if not output_file.exists():
            raise Exception("HH-blits failed. Protein is probably too long ?")
    return hhr_file


def call_trimal(input_file, output_file, trimal_gt, cons, colnumbering_file):
    print("calling trimal gt "+str(trimal_gt)+" cons "+str(cons)+" on "+str(input_file))
    fm.check_if_file_ok(input_file)
    cmd = "trimal -in "+str(input_file)+" -out "+str(output_file)+" -gt "+str(trimal_gt)+" -cons "+str(cons)+" -colnumbering > "+str(colnumbering_file)
    subprocess.Popen(cmd, shell=True).wait()


def call_reformat(input_file, output_file):
    print("will reformat "+str(input_file)+" thanks to Soeding's reformat.pl to "+str(output_file))
    fm.check_if_file_ok(input_file)
    call = "reformat.pl a3m fas "+str(input_file)+" "+str(output_file)+" -r"
    subprocess.Popen(call, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).wait()


def call_muscle(input_file, output_file=None):
    print("Calling Muscle")
    if output_file is None:
        output_file_name = '.'.join(str(input_file).split('.')[:-1])+"_muscle.fasta"
        output_file = pathlib.Path(output_file_name)
    muscle_cline = MuscleCommandline(input=str(input_file), out=str(output_file))
    print(muscle_cline)
    fm.check_if_file_ok(input_file)
    stdout, stderr = muscle_cline()
    return output_file


def call_muscle_profile(msa_file, seq_file, output_file):
    fm.check_if_file_ok(seq_file)
    muscle_cline = MuscleCommandline(profile=True, in1=str(msa_file), in2=str(seq_file), out=str(output_file), gapopen=-0.1)
    stdout, stderr = muscle_cline()
