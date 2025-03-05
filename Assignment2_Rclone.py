import subprocess
import shlex

class FileManager:
    def __init__(self, remote):
        self.remote = remote
        self.master_cloud = "dropbox:/Angela Huynh/8802/Assignment2"
        self.master_local = "/Users/ahuynh/Downloads/8802/TempData"

    def convert(self, file):
        '''
        - Creates the local and cloud file path equivalents of an input file
        Input: (str) single name of file
        Returns: (tuple) full paths of local file and cloud file
        '''
        local_filename = f"{self.master_local}/{file}"
        cloud_filename =  f"{self.master_cloud}/{file}"
        return local_filename, cloud_filename
    
    def uploadData(self, file):
        '''
        - Uploads a local file to the cloud
        Input:  (str) single name of local file

        rclone copy ~/downloads/4400_final "dropbox:/Angela Huynh/8802/4400" -P
        rclone copyto /Users/ahuynh/Downloads/8802/TempData/hello.py "dropbox:/Angela Huynh/8802/Assignment2/" -P
        '''
        res = self.convert(file) # this is converted cloud file
        local_filename, cloud_filename = res

        res = f'rclone copyto {local_filename} "{cloud_filename}" -P'

        res = shlex.split(res) # turning string into list, ignoring quoted content
        subprocess.run(res)
        print("Success! Uploaded {file} from {self.master_local} to {self.master_cloud}")


    def downloadData(self, file):
        '''
        Downloads cloud file to the local drive
        rclone copy "dropbox:/Angela Huynh/8802/AH_coverletter.pdf" ~/downloads/4400_final
        rclone copy google:MyDrive/Documents/document.pdf ~/Downloads

        '''
        res = self.convert(file)
        local_filename, cloud_filename = res
        res = f'rclone copy "{self.master_cloud}/{file}" {self.master_local}'
        print(res)
        res = shlex.split(res)
        subprocess.run(res)
        print("Success! Copied {file} from DropBox to {self.master_local}")



# Initializing objects
file_ob = FileManager("dropbox")
file_ob.uploadData("bye.py") # uploading bye.py from from the master_local directory
file_ob.downloadData("sam.py") # downloading sam.py from drop box directory to the master_local directory