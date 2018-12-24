#@Lonewolf v 0.1
import sys
import mmap
import re

#Uso: idMutaciones.py seqFile.fasta

programName = "idMutaciones"
programVersion = "0.1"
licenseVersion = " GNU v 3.0"
homepageProyect = "https://github.com/J0MS/idMutaciones_Tools"
programName.__repr__()
"Version".__repr__()
print(programName,"\n" ,"Version", programVersion,"\n" ,"license",licenseVersion,"\n" ,"Homepage:",homepageProyect)

try:
    if ((len(sys.argv) == 2)):

        fasta_File= sys.argv[1]
        #print("Args right",str(len(sys.argv)))
        sequence_Names = []
        sequences = []
        words = ''
        filepath = fasta_File
        with open(fasta_File, 'rb', 0) as file, mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as mapFile:
            line = mapFile.readline()
            line = line.decode('ISO-8859-1')
            cnt = 1
            while line:
                words += line
                match =re.findall(r'(>[A-Z]\w+)', line)
                if match:
                    sequence_Names.append(match)
                    #print(match)
                line = mapFile.readline()
                line = line.decode('ISO-8859-1')
                cnt += 1
            index= 0

        def find_between( s, first, last ):
            try:
                start = s.index( first ) + len( first )
                end = s.index( last, start )
                return s[start:end]
            except ValueError:
                return ""

        def find_between_r( s, first, last ):
            try:
                start = s.rindex( first ) + len( first )
                end = s.rindex( last, start )
                return s[start:end]
            except ValueError:
                return ""

        def compareSequences(seqName,a, b):
            index = 0
            pos = ""
            result = ""
            for x, y in zip(a, b):
                index += 1
                if x != y:
                    pos = index
                    gap,pos = isGap(a,b)
                    modification = checkMutationType(gap,a,b)
                    representation1=  x+str(index)+y
                    representation2= representation1.replace(" ", "")
                    wtCodon= codonIdexer(a,index)
                    mutantCodon = codonIdexer(b,index)

                    wtAmino = translate(wtCodon)
                    mutantAmino = translate(mutantCodon)
                    #mutant_is_synonymous= ""
                    mutant_is_synonymous= "sinónima," if translate(wtCodon)==translate(mutantCodon) else "no sinónima,"
                    aminos = wtAmino+mutantAmino
                    aminos = wtAmino+mutantAmino if '_' not in aminos else "(Corrimiento del ORF)"
                    mutantAminoRepresentation1= wtAmino+"-"+mutantAmino
                    mutantAminoRepresentation2=mutantAminoRepresentation1 if 'no sinónima' not in mutant_is_synonymous else wtAmino+str(index)+mutantAmino
                    efect= ""
                    if 'no sinónima' not in mutant_is_synonymous:
                       efect= "Nulo"
                    elif 'no sinónima' in mutant_is_synonymous:
                       efect= "Letal,"
                    else:
                       efect= "Leve"
                    print ("Mutation en",seqName,"posicion",index,"en las bases", x,y,"del tipo:",\
                    modification,"representacion:", representation2,",clase",mutant_is_synonymous,"codones(WT/MUT)",wtCodon,mutantCodon,\
                    "Amino(WT/MUT)",aminos,"representacion",mutantAminoRepresentation2)#"efecto:",efect)
            #return seqN,pos,baseX,baseY

        def codonIdexer(sequence,index):
            codons= []
            pointers = []
            pointer = 0
            key = ""
            for i in range(0, len(sequence), 3):
                temp = sequence[i:i + 3]
                codons.append(temp)
                pointers.append(i)
            dictPointers = {}

            for num in pointers:
                val = [num,num + 1,num +2]
                dictPointers[num] = val
            for k,v in dictPointers.items():
                for items in v:
                    if items == index:
                        pointer = int(k)
                        break

            indexKey = 0
            dictCodons = {}
            for x, y in zip(pointers,codons):
                dictCodons[x]=y

            for k,v in dictCodons.items():
                if k == pointer:
                    key = v
            return key

        def isGap(wtSeq,querySeq):
            gap = False
            index = 0

            if len(wtSeq) != len(querySeq):
                gap = True
                for x, y in zip(wtSeq,querySeq):
                    index += 1
                    if x != y:
                        break
            if index != 0:
                index = index
            return gap,index


        def checkMutationType(gap,wtSeq,querySeq):
            modificationClass=""
            if gap:
                modificationClass = "deleción"
            elif len(wtSeq) > len(querySeq):
                modificationClass = "adición"
            else:
                modificationClass = "Sustitución,"
            return modificationClass

        def translate(seq):

            table = {
                'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
                'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
            }
            protein =""
            if len(seq)%3 == 0:
                for i in range(0, len(seq), 3):
                    codon = seq[i:i + 3]
                    protein+= table[codon]
            return protein

        dictionary = {}
        words = words.replace('\r', '').replace('\n', '')
        index= 0
        for name in range(len(sequence_Names)-1):
            key  = ''.join(sequence_Names[index])
            str1 = ''.join(sequence_Names[index])
            str2 = ''.join(sequence_Names[index+1])
            values = find_between( words, str1,str2 )
            dictionary[key] = find_between( words, str1,str2 )
            index = index + 1

        #Agregada la cadena ">FIN" al archivo original para agregar al dicionario la ultima secuencia.
        dictionary[''.join(sequence_Names[len(sequence_Names)-1])] = find_between( words, ''.join(sequence_Names[len(sequence_Names)-1]),">FIN" )
        wildTypeSeq = dictionary.get(''.join(sequence_Names[0]))

        #Retiro ">FIN" al dicionario para evitar lecturas de secuencias con limites erroreas.
        dictionary.pop(">FIN")
        for k, v in dictionary.items():
            mutantSeq = dictionary.get(k)
            gap,pos = isGap(wildTypeSeq,mutantSeq)
            modification = checkMutationType(gap,wildTypeSeq,mutantSeq)
        #    sN,i,x,y = compareSequences(k,wildTypeSeq,mutantSeq)
            if gap:
                representation1=  "∆" + mutantSeq[pos]+str(pos)
                representation2= representation1.replace(" ", "")
                print("*****Gap en la secuencia",k,"en la posicion",pos,"****************","tipo de mutacion,",modification,",  representacion:", representation2 )
            else:
                #sN,i,x,y = compareSequences(k,wildTypeSeq,mutantSeq)
                #print ("Mutation en la secuencia",sN,"en la posicion",i,"en las bases", x,y )
                compareSequences(k,wildTypeSeq,mutantSeq)
                #print(tmp, modification)
    else:
        print ("No arguments detected","\n","Usage:","idMutaciones.py seqFile.fasta")


except:
    print ("Faltal, I/O Failed!")
    print("Args failed",str(len(sys.argv)))
