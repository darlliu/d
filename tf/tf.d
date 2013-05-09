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
    double up, down;
    double pval;
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
                tfmap [(*temp).TFname] ~= temp;
                probemap [(*temp).probeid] ~= temp;
            }
       }
    };
    void gen_avg()
    {
        //Generate the average behavior a central TF
    };
    void gen_tfs()
    {
        //iterate through the tfs and do the fisher test
    };
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
    foreach (string key, datapt*[] _val; mymap.tfmap)
    {
        io.writefln("The key is this: %s",key);
        foreach (datapt* _val2; _val)
        {
            auto val = *_val2;
            io.writefln ("The three values again are %s, %s, %s",
                    val.TFname,val.probeid, val.genesym);
            io.writefln ("The numbers are %s, %s, %s", 
                    val.wt_mean, val.exp_mean, val.pval);
            //assert (key == val.TFname);
        }
    }
}
