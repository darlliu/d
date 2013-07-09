/*
 *Acquire the following:
 *1, given one table of Cyber-T (probeid, genesym, wt-mean, exp-mean, pval)
 *2, given another one with (TF, genesym, probeid)
 *
 *Construct the following:
 *1, a struct with the info of the above, the genesym from Cyber-T is underscored
 *2, a two way hash array for lookup both in terms of TF and in terms of probeid
 *
 *Perform the following:
 *1, Fisher exact test on UP and DOWN p values compared with the average of all TFs.
 *2, t test again against the above
 *
 *Output the following:
 *1, (TFname, pval_up, pval_down, expression (if can be found))
 *2, permute the pval and target
 *3, also relative numbers
 *
 *for top 100, and all of the TFs.
 */

import io=std.stdio;
import reg= std.regex;
import str=std.string;
import std.conv;
struct datapt {
    string TFname, probeid, genesym, _genesym;
    double wt_mean, exp_mean, pval;
}
struct tfpt {
    string name;
    uint up, down, insig,targets;
    real pval_up, pval_down;
    double prec_up, prec_down;
}

static pure real c (in uint n, in uint c) 
{
    //produces cCn in a relatively safe way
    //in general n>>c, no need to check.
    if (n==0||c==0) return 1;
    real acc= 1;
    for (ulong k=n, m=c; k>=n-c+1 || m>=1; k--, m--)
    {
        if(k>=n-c+1) acc*=k;
        if(m>=1) acc/=m;
        //io.writefln("Here we got k, m, acc as %f, %f, %f", k, m, acc );
    }
    return acc;
}
static pure void fisher(in tfpt wt, ref tfpt exp)
{
    //now get up
    exp.pval_up= c(wt.up+exp.up,wt.up)*(
            c(wt.down+wt.insig+exp.down+exp.insig,wt.down+wt.insig)
            /
            c(wt.up+wt.down+wt.insig+exp.up+exp.down+exp.insig,wt.up+wt.down+wt.insig)
                );
    //io.writeln(to!string(exp.pval_up));
    //now get down
    exp.pval_down= c(wt.down+exp.down, wt.down)*(
            c(wt.up+wt.insig+exp.up+exp.insig,wt.up+wt.insig)
            /
            c(wt.up+wt.down+wt.insig+exp.up+exp.down+exp.insig,wt.up+wt.down+wt.insig)
                );
    //io.writeln(to!string(exp.pval_down));
}
static string join (T) (in T[] l, in string set="\t")
{
    auto acc="";
    foreach (it; l)
    {
        acc~=to!string(it)~"\t";
    }
    return acc[0 .. $-1];
}
class tfmaps
{
    this ()
    {
    };
    this (in string tfname,in string expname)
    {
        this.bootstrap(tfname);
        this.load(expname);
    };
    void load (in string expname)
    {
        auto f = io.File(expname, "r");
        scope(failure){ io.writeln("Error loading file!!"); };
        scope(success){ io.writeln("Finished loading the expression data"); };
        foreach (string line; io.lines(f))
        {
            auto words = reg.split(line,reg.regex("\t"));
            string probe = str.strip(words[0]);
            //io.writefln("The words are %s", words);
            if (probe !in probemap) {
                //io.writefln("Probe %s was not in the map!", probe);
                continue;
            }
            auto _temp = probemap[probe];
            foreach (datapt* temp; _temp)
            {
                (*temp)._genesym = str.strip(words[1]);
                //io.writefln("Parsed sym is %s", (*temp)._genesym);
                (*temp).wt_mean =  to!double(str.strip(words[2]));
                //io.writefln("Parsed wtmean is %s", (*temp).wt_mean);
                (*temp).exp_mean =  to!double(str.strip(words[3]));
                //io.writefln("Parsed expmean is %s", (*temp).exp_mean);
                (*temp).pval = to!double(str.strip(words[4]));
                //io.writefln("Parsed pval  is %s", (*temp).pval);
                tfmap [(*temp).TFname] ~= temp;
            }
        }
    };
    void bootstrap(in string tfname)
    {
        auto f = io.File(tfname,"r");
        scope(failure){ io.writeln("Error opening file!!"); };
        scope(success){ io.writeln("Finished bootstraping the datapoints!"); };
        foreach ( string line;  io.lines(f) )
        {
            auto words = reg.split(line, reg.regex(","));
            if (words.length < 2) throw new Exception("Read Opening");
            auto _words = reg.split(words[0],reg.regex("\t"));
            //io.writeln("Words and _words are ",words, _words);
            auto temp = new datapt;
            //with (temp)
            {
                (*temp).TFname = str.strip(_words[0]);
                (*temp).genesym = str.strip(_words[1]);
                (*temp).probeid = str.strip(words[1]);
                //io.writefln("Here we got %s, %s, %s",
                        //(*temp).TFname , (*temp).genesym, (*temp).probeid);
                probemap [(*temp).probeid] ~= temp;
            }
       }
    };
    void gen_tfs(in double CUTOFF=5e-2)
    {
        //iterate through the tfs and do the fisher test
        uint up_acc , down_acc , insig_acc ;
        up_acc = down_acc = insig_acc = 0;
        foreach (it; tfmap)
        {
            tfpt temp;
            with (temp)
            {
                name= (*it[0]).TFname;
                targets=to!uint(it.length);
                up = down = insig = 0;
                //io.writefln("At %s, has this many targets: %d", name, it.length);
                foreach(itt; it)
                {
                    if ((*itt).pval < CUTOFF )
                        (*itt).wt_mean>(*itt).exp_mean?down++:up++;
                    else insig++;
                }
                prec_up=to!double(up)/(up+down+insig);
                prec_down=to!double(down)/(up+down+insig);
                up_acc+=up;
                down_acc+=down;
                insig_acc+=insig;
                //io.writefln("At %s, up %d, down %d, insig %d", name, up, down, insig);
            }
            tfs~=temp;
        }
        with (AVG)
        {
            name="AVERAGE_TF";
            up=up_acc/to!uint(tfmap.length);
            down=down_acc/to!uint(tfmap.length);
            insig=insig_acc/to!uint(tfmap.length);
            io.writefln("Length %d, insig %d", tfmap.length, insig);
        }
    };
    void fisher_test()
    {
            //io.writefln("Now doing fisher");
        foreach (ref it; tfs)
        {
            fisher(AVG,it);
            //io.writefln("Changed pvals to %f, %f", it.pval_up, it.pval_down);
        }
    };
    void to_print(in string ofname="output.tsv")
    {
        auto f=io.File(ofname,"w");
        f.write(join([
                    "TFname", "Targets_Up","Targets_Down", "Targets_Neither","Targets_total",
                    "pVal_Up", "pVal_Down","%Up","%Down","\n"
                    ]));
        foreach(it; tfs)
        {
            with(it)
            {
                f.write(name~"\t"~join([
                        up, down,insig,targets,pval_up,pval_down, prec_up,prec_down
                        ])~"\n");
            }
        }
        f.close();

    }
    datapt*[][string] probemap;
    datapt*[][string] tfmap;
    tfpt AVG;
    tfpt[] tfs;
}
int main(string[] args)
{
    auto mymap = new tfmaps;
    if (args.length < 3) return -1;
    try
    {
        mymap.bootstrap(args[1]);
        mymap.load(args[2]);
        mymap.gen_tfs();
        mymap.fisher_test();
        mymap.to_print();
    }
    catch
    {
        return 1;
    }
    return 0;
}

unittest
{
    auto mymap = new tfmaps;
    mymap.bootstrap("./tf_affy_map_hg.txt");
    mymap.load("./exported.txt");
    //foreach (string key, datapt*[] _val; mymap.tfmap)
    //{
        ////io.writefln("The key is this: %s",key);
        //foreach (datapt* _val2; _val)
        //{
            //auto val = *_val2;
            ////io.writefln ("The three values again are %s, %s, %s",
                    ////val.TFname,val.probeid, val.genesym);
            ////io.writefln ("The numbers are %s, %s, %s", 
                    ////val.wt_mean, val.exp_mean, val.pval);
            //assert (key == val.TFname);
            //assert (val._genesym!="");
        //}
    //}
    mymap.gen_tfs();
    mymap.fisher_test();
    io.writefln("The average have stats %d, %d, %d"
            ,mymap.AVG.up, mymap.AVG.down, mymap.AVG.insig);
    foreach(it; mymap.tfs)
    {

        io.writefln ("Genesym %s:%d targets, %d up, %d down,%f pre up %f pre down,  %f pvalup %f pval down",
                it.name, it.targets, it.up,it.down,it.prec_up, it.prec_down,it.pval_up,it.pval_down
                );
    }
    io.writefln("C5,2 is %f", c(5,2));
    assert( to!int(c(5,2))==10);
    io.writefln("C6,3 is %f", c(6,3));
    assert( to!int(c(6,3))==20);
    io.writefln("C8,2 is %f", c(8,2));
    assert( to!int(c(8,2))==28);
    io.writefln("join is this %s", join([1,2,3]));
    mymap.to_print();
    //assert( ft(3)==6);
    //assert( ft(10)==3628800);
}
