# Steps to install Rceattle from the source file

# Step 1 - Download Rceattle_0.0.0.9000.tar.gz from:
# https://drive.google.com/drive/folders/1sHe_KxvKZi7UyjnWB4Oz9HmKY6CWGWoM

# Step 2 - set "file_directory" to the directory where the download is
file_directory <- "your_download_directory"
path_to_file <- paste(file_directory, "/", "Rceattle_0.0.0.9000.tar.gz", sep = "")

# Step 3 - Install
install.packages(path_to_file, repos = NULL, type="source")
