from sys import argv    
import os.path
    
if __name__ == '__main__':
	mappingfl = argv[1]
	
	filenms = []
	filenms.append(mappingfl)
	for i in range(2, len(argv)):
		filenms.append(argv[i])
		
	open("dataFilesforJS.txt",'w').writelines([f+'\n' for f in filenms])