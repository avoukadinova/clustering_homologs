
## A) Overview

When given a plethora of gene sequences (to be exact, over 100,000) from various organisms, the mere idea of determining which sequences are most similar can make one’s head (and computer) crash. Comparative biology and the ability to recognize homologous (similar in position, structure, and evolutionary origin but not necessarily in function) genes (i.e., homologs) is crucial to genome annotation, protein structure, and phylogenetics. This large set of next-generation sequencing data requires high-speed, efficient, and convenient bioinformatics algorithms to accurately execute clustering sequences. 

The current standard to clustering homologs is USEARCH due to its impecable run time. Although faster, USEARCH often fails to cluster sequences that are actually homologous due to the nature of the algorithm. Our solution to USEARCH is an executable algorithm called clust_horo that retains the desired “fast” speed of USEARCH, but also removes the bias and iteratively performs clustering on the UCLUST output to produce a more accurate homolog clustering of the data, thus removing any unnecessary clusters. Effectively, clust_horo is still able to handle large data sets of 80,000 to 200,000 gene sequences while maintaining a relatively low run time.

For the clust_horo algorithm to work, Biopython and USEARCH must be installed on your machine. A path to USEARCH must be added to your local bin so that USEARCH can be run from anywhere. Please follow ALL of the installation instructions below to ensure clust_horo works with no errors or issues. It should be noted that clust_horo is an executable, so all work is done through the command line. The instructions below will not work from the Windows Command Prompt.

## B) Installing Biopython

For Biopython to be installed, you must first have Python installed. As there are many ways to do this, we leave it up to you to figure out. 
  
```{r eval=FALSE,echo=TRUE}
pip install biopython
```

If this doesn't work for you, other options for installing Biopython can be found here: <http://biopython.org/DIST/docs/install/Installation.html>

## C) Installing USEARCH

### 1) Download USEARCH

Please visit <https://www.drive5.com/usearch/download.html> and choose your specific operating system. From there, you will get an email from Robert Edgar with a link. Click on this link. USEARCH will download automatically. The file that is downloaded is an executable, so you cannot open it or click on it!

### 2) Move the executable to your desired directory

Let's say the file I just downloaded is in my downloads folder (for the sake of these instructions, let's say the path is /Users/Name/Downloads/). Move the executable to your desired location where you want it to permanently be. This can be done by clicking and dragging or through the command line.
Let's say I want to move this file to a folder named USEARCH (path: /Users/Name/Desktop/USEARCH/). Run the following command.
  
  
```{r eval=FALSE,echo=TRUE}
mv /Users/Name/Downloads/usearch11.0.667_i86osx32 /Users/Name/Desktop/USEARCH
```
  
### 3) Go into the USEARCH directory  

From here on out, everything is done through the command line. Go into the directory where you just moved your executable to.    
  
  
```{r eval=FALSE,echo=TRUE}
cd /Users/Name/Desktop/USEARCH
```

### 4)  

Change the access permission of the executable.  
  
  
```{r eval=FALSE,echo=TRUE}
chmod 744 usearch11.0.667_i86osx32
```

### 5)  

Create a path to USEARCH in your local bin. You will be prompted to input your password.   
   
  
```{r eval=FALSE,echo=TRUE}
sudo ln -s /Users/Name/Desktop/USEARCH/usearch11.0.667_i86osx32 /usr/local/bin/usearch/
```

You can now run USEARCH from anywhere on your machine!! If this did not work, USEARCH installation instructions can be found here: <https://www.drive5.com/usearch/>. An example of a USEARCH command you can run is

```{r eval=FALSE,echo=TRUE}
usearch -cluster_fast input_seq.fasta -id 0.8 -centroids centroid.txt -clusters cluster.txt
``` 

USEARCH documentation can be found here: <https://www.drive5.com/usearch/manual/>

## D) Downloading clust_horo  
  
The clust_horo algorithm can be found here: <https://github.com/avoukadinova/clustering_homologs/>

## E) Running clust_horo

Once you have downloaded clust_horo, go into the directory where the executable is. To run it, one of two commands should work.

```{r eval=FALSE,echo=TRUE}
./clust_horo.py -i inputseq.fasta -t 0.6 -m 10
``` 

or

```{r eval=FALSE,echo=TRUE}
python3 clust_horo.py -i inputseq.fasta -t 0.6 -m 10
``` 
For a list of all the flag parameters, run the following command. 

```{r eval=FALSE,echo=TRUE}
./clust_horo.py -h
```

## F) Reclustering

Once you run clust_horo, the only output will be the initial clusters made by USEARCH, the centroids text file, and a homo_out.txt file. You, as the user, have two options:

1) Recluster manually
2) Run reclustering.py

To view the homologous clusters, open the homo_out.txt file. To automatically recluster, download recluster.py, add it to the same directory as clust_horo.py, and run the following command. 

```{r eval=FALSE,echo=TRUE}
./recluster.py
```
or

```{r eval=FALSE,echo=TRUE}
python3 recluster.py
```

recluster.py uses the homo_out.txt file to recluster. The output is an updated centroids text file and the accurate clusters in FASTA format.

