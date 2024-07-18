# Uploading data to GEO
Author: Mary Allen and Lynn Sanford

When publishing on NIH grants, Before you publish you have to upload data to the Gene Expression Omnibus (GEO). GEO [https://www.ncbi.nlm.nih.gov/geo/](https://www.ncbi.nlm.nih.gov/geo/) contains all the project metadata, and it links to the raw fastq files that are stored in the Short Read Archive (SRA) database.

## Log into GEO
[Sign up or log into GEO](https://account.ncbi.nlm.nih.gov/signup/?back_url=https://www.ncbi.nlm.nih.gov/geo/submitter/) to get a personal folder where your data will go.
    - I’d suggest using your eRA commons ID to log in if you have one. It's assigned by your university and is used by your PI for various reasons. 

Once you are logged in, follow [these directions](https://www.ncbi.nlm.nih.gov/geo/info/submissionftp.html). This should be done before anything else because GEO will provide you an upload folder and password.
    - For example, mine is: `uploads/maallen3_KPEXeeOm`

## Follow submission steps
Detailed instructions for uploading are provided [here](https://www.ncbi.nlm.nih.gov/geo/info/seq.html).

### Fill out metadata spreadsheet
The instructions above link to the template spreasheet to use for metadata. It is extremely detailed, and for good reason. The more detailed metadata you can provide about your data, the more confident others will be in using it in the future.

### Create your submission folder
This contains all the raw and processed files you'd like to upload.

### Upload data
Go to your submission folder on your server/supercomputer, then follow these steps:
- Type `sftp geoftp@sftp-private.ncbi.nlm.nih.gov`
- When prompted, `password: <geo password they gave you>`
- You are now on a ftp server. It looks the same as a bash server but doesn’t use all the same commands.
- `cd uploads/<username_randomprojectcodetheygiveyou>`
    - In my case, `cd uploads/maallen3_KPEXeeO` 
- For one file:
    ```
    mkdir new_geo_submission
    cd new_geo_submission
    put <file_name>
    ```
- For one directory full of files
    ```
    mkdir mydir
    put -r mydir 
    ```
    - This command is grabbing the data from your current working directory on your server/supercomputer and putting it on the ftp server

## Errors
If you get the message "The RSA host key for sftp-private.ncbi.nlm.nih.gov has changed, and the key for the corresponding IP address", it means you need to delete the previous line in your known_hosts file, which was created the last time you logged into the NIH FTP server.

```
vim /Users/<username>/.ssh/known_hosts
```
Find and delete the line with `sftp-private.ncbi.nlm.nih.gov`.
- Remember `dd` in vim deletes the whole line (in command mode)



