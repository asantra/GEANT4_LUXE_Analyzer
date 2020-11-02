### script to submit condor jobs
### usage:
#### python sendCondorJobs.py -l <fileList> -- short (for short job)
import os, sys, time
import argparse

def main():
    topDir = "/afs/desy.de/user/s/santraar/private/Geant4_File_Analyzer/BackgroundAnalyzer/"
    parser = argparse.ArgumentParser()
    parser.add_argument('--short', dest='shortJob', action='store_true')
    parser.add_argument('-l', dest='fileList', action='store', type=str, help="The file list")
    args = parser.parse_args()
    
    inputFileList = args.fileList
    needShort     = args.shortJob

    numFile = 0
    with open(inputFileList) as inFile:
        for lines in inFile.readlines():
            rootFileName = lines.rstrip()
            ### bash file writing
            fileSH    = open('runBKGMaster.sh')
            fileOutSH = open('runBKG_v'+str(numFile)+'.sh', 'w')
            for line in fileSH.readlines():
                if 'RRRR' in line:
                    line = line.replace('RRRR', rootFileName)
                if 'VVVV' in line:
                    line = line.replace('VVVV', str(numFile))
                fileOutSH.write(line)
            fileOutSH.close()
            
            ### prepare the jdl script
            fileJdl    = open('condorJobsMASTER.jdl')
            fileOutJdl = open('condorJobs_v'+str(numFile)+'.jdl', 'w')
            for line in fileJdl.readlines():
                if 'XXX' in line: line = line.replace('XXX', 'runBKG_v'+str(numFile)+'.sh')
                #### not needed as everything in cvmfs ###
                if 'ZZZ' in line: line = line.replace('ZZZ',topDir+"process_track_tree_draw_v4.C")
                if 'MMC' in line: line = line.replace('MMC',topDir+"MHists.C")
                if 'MMH' in line: line = line.replace('MMH',topDir+"MHists.h")
                if 'AAA' in line: line = line.replace('AAA','stdout_v'+str(numFile)+'.out')
                if 'BBB' in line: line = line.replace('BBB','stderr_v'+str(numFile)+'.err')
                if 'CCC' in line: line = line.replace('CCC','batchlog_v'+str(numFile)+'.log')
                if 'HHH' in line: line = line.replace('HHH','allBackgroundFiles_v'+str(numFile)+".root")
                    
                fileOutJdl.write(line)
            fileOutJdl.close()
            
            ### submit the jobs
            os.system("condor_submit condorJobs_v"+str(numFile)+".jdl")
            
            numFile += 1
            if(needShort and numFile>1):
                break
            
            
if __name__=="__main__":
    start = time.time()
    main()
    print "The time taken: ", time.time() - start, " s"
