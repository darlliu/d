#!/auto/igb-libs/linux/centos/6.x/x86_64/pkgs/python/2.6.8/bin/python
f= open("./exported_sirt1_ko.txt","r")
g= open("./exported_sirt1_ko_2.txt","w")
for line in f:
    words=line.split("\t")
    words[0]=words[0].strip("\"").upper()
    g.write("\t".join(words));
g.close()
f.close()

f= open("./exported_sirt1_oe.txt","r")
g= open("./exported_sirt1_oe_2.txt","w")
for line in f:
    words=line.split("\t")
    words[0]=words[0].strip("\"").upper()
    g.write("\t".join(words));
g.close()
f.close()
