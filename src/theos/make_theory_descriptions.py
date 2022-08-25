import glob

list=[]
theos=glob.glob('../theos/th_*.[Ff]*')
theos+=glob.glob('../unused_theos/th_*.[Ff]*')

theos.sort()    
desc=[]
desclist=[]
for theo in theos:
    if theo.strip()[-1]=='~':
        continue
    file = open(theo)
    lines=file.readlines()
    file.close()
    parnam=[]
    func=-1
    for l in range(len(lines)):
        if lines[l].lstrip().lower().startswith('function') and func==-1: # only the first is the func definition
            name=lines[l].lstrip().lower().split()[1].split('(')[0]
            func=l
        if lines[l].lstrip().lower().startswith('parnam'):
            parnam.append(lines[l].lstrip().lower().strip())
    #anchor in wiki
    #<div id="NameOfAnchorHere">optional text</div>
    desc.append('<div id="'+name+'">'+"'''"+name+"'''"+'</div>'+'    filename '+theo)
    desclist.append('[[#'+name+']]')
    while (lines[func+1].lstrip()+' ')[0]=='!':
        
        desc.append(lines[func+1].strip()[1:].replace('!','',3).replace('===','').replace('---','').strip())
        func+=1
    desc+=parnam
    desc.append('----')
f=open('theory_description.txt','w')
f.writelines(['#'+line+'<br />\n'for line in desclist])
f.writelines([line+'<br />\n'for line in desc])
f.close()

