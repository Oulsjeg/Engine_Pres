
## **Guide to Interact with Git Branches**
=====================
### **In Github Desktop**
=====================\
Please read this:\
https://help.github.com/en/desktop/contributing-to-projects/making-changes-in-a-branch


=====================
### **In gitkraken**
=====================\
Please read this:\
https://support.gitkraken.com/working-with-repositories/branching-and-merging/

====================
### **In terminal**
====================

Create branch:\
`$git branch <branch_name>`

Routine fetch and checkout your branch #"checkout" means that you switch to current working status to the branch <branch_name>:\
`$git fetch && git checkout <branch_name>`

In main branch, list all the branches names:\
`$git branch -a`

Stage all files and Commit to a branch, use git add <file_name> if you  just want to stage specific file:\
`$git checkout <branch_name>`\
`$git pull`\
`$git add .`\
`$git commit -m "input whatever commit message"`\
`$git push origin <branch_name>`

Switch back to the master branch:\
`$git checkout master`

Before you switch to your master branch:\
Note that if your working directory or staging area has uncommitted changes that conflict with the branch you’re checking out, Git won’t let you switch branches.\
Merge a branch with the master branch:\
`$git checkout master`\
`$git merge <branch_name>`

Delete a branch\
`$git branch -d <branch_name>`\
`$git branch -D <branch_name>`


=====================
  ### **More tricks**
=====================

Combine $git branch <branch_name> && git checkout <branch_name> by one line:\

`$git checkout -b <branch_name>`

Combine $git add . && git commit -m 'messages'

`$git commit -am 'messages'`

What if there is a merge conflict\
E.g., You changed a=[2 3 4 5] to a=[1 3 4 5] at the master branch and someone else changed a=[2 3 4 5] to a=[2 3 4 6] at a <branch_name> branch. Then you guys merge. What will happen. Anything that has merge conflicts and hasn’t been resolved is listed as unmerged. Git adds standard conflict-resolution markers to the files that have conflicts, so you can open them manually and resolve those conflicts. Your file contains a section that looks something look like this:
>>> `<<<<<<<` HEAD:index.py\
>line 30: a=[1 3 4 5]\
>=======\
>line 32:\
>a=[2 3 4 6]\
>>> `>>>>>>>` branch_name:index.py

If you want to use a graphical tool to resolve these issues, you can run git mergetool, which fires up an appropriate visual merge tool and walks you through the conflicts:\
`$git mergetool`\
then you runthis line to verify that all conflicts have been resolved:\
`$git status`
>>>On branch master\
All conflicts fixed but you are still merging.\
  (use "git commit" to conclude merge)\
  Changes to be committed: modified:   index.py'

More reading about resolving merge conlicts:\
https://www.atlassian.com/git/tutorials/using-branches/merge-conflicts
https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging
