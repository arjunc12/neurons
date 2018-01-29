import os

ASC_DIR = 'sandbox'

def write_swc(lines, outname):
    f = open(outname, 'w')
    for line in lines:
        outstr = ' '.join(line)
        f.write('%s\n' % outstr)
    f.close()

def split_swc(fname):
    curr_lines = None
    curr_file = 0
    assert fname[-4:] == '.swc'
    prefix = fname[:-4]
    for line in open(fname):
        if line[0] == '#':
            continue
        else:
            line = line.strip('\n')
            line = line.split()
            assert len(line) == 7
            if line[-1] == '-1':
                if curr_lines != None:
                    outname = prefix + str(curr_file) + '.swc'
                    write_swc(curr_lines, outname)
                curr_lines = [line]
                curr_file += 1
            else:
                curr_lines.append(line)

def convert_asc():
    for fname in os.listdir('.'):
        prefix = fname[:-4]
        suffix = fname[-4:]
        print prefix, suffix
        if suffix == '.ASC':
            os.system('./NLMorphologyConverter %s.ASC %s.swc SWC' % (prefix, prefix))
            split_swc(prefix + '.swc')

def main():
    os.chdir('boutons/asc_files')
    #os.system('./NLMorphologyConverter G419g#1-E.ASC G419g#1-E.swc SWC')
    #split_swc('G419g#1-E.swc')
    convert_asc()

if __name__ == '__main__':
    main()
