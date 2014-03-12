import os
import sys



def getLibs(directory):
    libdir = {}
    for filename in os.listdir(directory):
        id = filename.partition('_')[0]
        libdir[id] = filename
    return libdir


def list():
    output = os.path.join(
        os.path.dirname(
            os.path.dirname(
                os.path.abspath(__file__))),
        'RNASeq.log')
    directory = '/home/vodkinlab/RNASeq/raw_sequences/'
    print directory
    libraries = getLibs(directory)
    i = 1
    ofh = open(output, 'w')
    while True:
        id = "R%s" % str(i).zfill(2)
        print id
        try:
            print libraries[id]
            ofh.write(libraries[id] + '\n')
            i += 1
        except:
            break
    ofh.close()
    return


def main():
    list()

if __name__ == '__main__':
    main()
