---
title: "Installing clust_horo"
output:
  pdf_document: default
  html_document: default
---

## A) Overview

For the clust_horo algorithm to work, Biopython and USEARCH must be installed on your machine. A path to USEARCH must be added to your local bin so that USEARCH can be run from anywhere on your machine. Please follow ALL of the installation instructions below to ensure clust_horo works with no errors or issues. It should be noted that clust_horo is an executable, so all work is done through the command line. The instructions below will not work from the Windows Command Prompt.

## B) Installing Biopython

For Biopython to be installed, you must first have Python installed. As there are many ways to do this, we leave it up to you to figure out. 

### Linux and MAC  
  
```{r eval=FALSE,echo=TRUE}
pip install biopython
```

If this doesn't work for you, you can try other options for installing Biopython here: <http://biopython.org/DIST/docs/install/Installation.html>

## C) Installing USEARCH

### 1) Download USEARCH

Please visit <https://www.drive5.com/usearch/download.html> and choose to download USEARCH for your specific operating system. From there, you will get an email from Robert Edgar with a link. Click on this link. USEARCH will download automatically. The file that is downloaded is an executable, so you cannot open it or click on it!

### 2) Move the executable to your desired directory

Let's say the file I just downloaded is in my downloads folder (for the sake of these instructions, let's say the path is /Users/Name/Downloads/). Move the executable to your desired location where you want it to permanently be. You can do this by clicking and dragging or through the command line.
Let's say I want to move this file to a folder named USEARCH (with this path: /Users/Name/Desktop/USEARCH/). The following command will do that.  
  
  
```{r eval=FALSE,echo=TRUE}
mv /Users/Name/Downloads/usearch11.0.667_i86osx32 /Users/Name/Desktop/USEARCH
```
  
### 3) Go into the USEARCH directory  

From here on out, everything is through the command line. Go into the directory where you just moved your executable to.    
  
  
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

You can now run USEARCH from anywhere on your machine!! If this did not work, USEARCH installation instructions can be found here: <https://www.drive5.com/usearch/>
  
## D) Downloading clust_horo  
  
The clust_horo algorithm can be found here: <https://github.com/avoukadinova/clustering_homologs/>