--------------------------------------

To create a new repository

ssh -i ~/.ssh/my-git-key eng-milab@git.csx.cam.ac.uk create medimag/NEWREPONAME

mkdir NEWREPONAME
cd NEWREPONAME
git init
# this just makes sure there's at least one file ...
touch README  
git add README
git commit -m 'initial commit' -a
git remote add origin eng-milab@git.csx.cam.ac.uk:medimag/NEWREPONAME
git push origin master

---------------------------------------

To update repository

git add . # when in dir where you want to add all files or specify by name
git commit -m 'message' -a
git push origin master

-----------------------------------------

To update to a previous repository

git log path/to/file # work out which revisions to revert too
git checkout <commit> path/to/file # revert - <commit> is the descriptor tag


-----------------------------------------

Update a file to HEAD

git checkout path/to/file ~ will overite any uncommited local changes

-----------------------------------------

To clone the repository

git clone eng-milab@git.csx.cam.ac.uk:medimag/rosePhD.git NEWREPOFOLDERNAME

-----------------------------------------

To merge a detached head

git commit -m "commit the detached head"
git checkout -b BRANCH_NAME
git checkout master
git merge merge BRANCH_NAME # conficts result
git status # not the confilted files
git checkout BRANCH_NAME file_name1 file_name2 ...
git add file_name1 file_name2 ...
git commit -m "another commit message" -a
git push origin master

-------------------------------------------

To see changes in files

git diff HEAD optionalFileName.cxx # diff between current dir (uncommitted) and pervious commit
git diff HEAD^ # diff between current dir and two commits ago
git diff HEAD^ HEAD # diff between two commits. 

----------------------------------------

Git - branch the return to origional

git branch BRANCH_NAME
git commit . -m 'message'
git 


————————————————————

Split a git Repository

1. clone the repo - then filter to keep only one directory
git clone ssh://eng-milab@git.csx.cam.ac.uk/medimag/rosePhD.git freshrepo
git filter-branch --prune-empty --subdirectory-filter FOLDER_TO_KEEP BRANCH_TO_KEEP

2. make a new repo 
ssh -i ~/.ssh/my-git-key eng-milab@git.csx.cam.ac.uk create medimag/NEWREPONAME
cd NEWREPONAME
git init
git push origin master

3.change the filtered repo to point to the new repo
git remote set-url origin eng-milab@git.csx.cam.ac.uk:medimag/NEWREPONAME
git push origin master


