Used to download Argo data in a manner that allows for SCV processing
The original Argo data as downloaded from the US GODAE FTP is archived under the folder "/pub/":

/pub/outgoing/argo/geo/pacific_ocean/
/pub/outgoing/argo/geo/atlantic_ocean/
/pub/outgoing/argo/geo/indian_ocean/

----------------------------

Downloading Argo data from FTP server "usgodae.org/pub/outgoing/argo/"

Commands:

go to argo data directory:

i.e. for our systems...
cd /data/project1/data/argo

download all basins together:

wget -cNr -nH ftp://www.usgodae.org/pub/outgoing/argo/geo

or download separately by basin:
         
wget -cNr -nH ftp://www.usgodae.org/pub/outgoing/argo/geo/pacific_ocean
wget -cNr -nH ftp://www.usgodae.org/pub/outgoing/argo/geo/atlantic_ocean
wget -cNr -nH ftp://www.usgodae.org/pub/outgoing/argo/geo/indian_ocean

----------------------------
      
NOTE: the download preserves the original file structure:
./pub/outgoing/argo/geo/pacific_ocean/...
etc., which allows *updating* argo files *without* downloading the full dataset again, using the original time-stamp of the files
    
---------------------------- 
 
For final submission of Argo SCVs...   
Last Download: February 11 2020.
