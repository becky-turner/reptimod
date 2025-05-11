# Create github feature branches
#install.packages("gert")
library(gert)
# Create new branch
git_branch_create("dev", checkout = TRUE)

# Check we're on the branch
git_push()

# Push it up to GitHub
git_push(set_upstream = TRUE) # TRUE resets the automatic new push location to the new branch.
