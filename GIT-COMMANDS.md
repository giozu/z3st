<!-- # --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- -->

# Git Commands Cheat Sheet

A compact collection of useful Git commands for daily development, collaboration, and repository maintenance.

---

## Configuration

```bash
# Set user name and email (global)
git config --global user.name "Name"
git config --global user.email "your@email.com"

# Enable colored output
git config --global color.ui auto

# List all current configurations
git config --list

# Get a specific config value
git config --get <key>

# Initialize a new Git repository
git init

# Clone a remote repository
git clone <repository_url>

# Clone a specific branch
git clone --branch <branch_name> <repository_url>

# Check file status (tracked, untracked, staged)
git status

# Add a specific file to the staging area
git add <file>

# Add all modified files
git add .

# Commit staged changes with a message
git commit -m "Commit message"

# Amend the last commit message
git commit --amend -m "New message"

# Show unstaged differences
git diff

# Show staged differences
git diff --staged

# List local branches
git branch

# List local and remote branches
git branch -a

# Create a new branch
git branch <branch_name>

# Switch to an existing branch
git checkout <branch_name>

# Create and switch to a new branch
git checkout -b <branch_name>

# Show all local branches with their tracking info
git branch -vv

# Publish the new branch
git push -u origin <branch_name>

# Rename current branch
git branch -m <new_name>

# Delete a local branch (safe)
git branch -d <branch_name>

# Force delete a local branch
git branch -D <branch_name>

# Merge a branch into the current one
git merge <branch_name>

# Rebase current branch on another
git rebase <base_branch>

# Interactive rebase (to edit/reorder commits)
git rebase -i <base_branch>

# Show all remotes
git remote -v

# Add a new remote
git remote add <name> <url>

# Remove a remote
git remote rm <name>

# Rename a remote
git remote rename <old> <new>

# Fetch updates from remote (no merge)
git fetch

# Pull (fetch + merge)
git pull <remote> <branch>

# Pull with rebase
git pull --rebase <remote> <branch>

# Push local changes to remote
git push <remote> <branch>

# Delete a remote branch
git push <remote> --delete <branch>

# (alternative)
git push <remote> :<branch>

# Stash uncommitted changes
git stash

# Apply and remove the last stash
git stash pop

# Apply without removing from stash list
git stash apply

# List stashed changes
git stash list

# Drop a specific stash
git stash drop <stash_id>

# Undo local modifications (reset to HEAD)
git reset --hard HEAD

# Undo changes in a specific file
git checkout -- <file>

# Remove untracked files
git clean -f

# Remove untracked files and directories
git clean -fd

# Create a lightweight tag
git tag <tag_name>

# Create an annotated tag
git tag -a <tag_name> -m "Message"

# Push a specific tag
git push <remote> <tag_name>

# Push all tags
git push <remote> --tags

# Prune unreachable objects (cleanup dangling commits/objects)
git prune

# Fetch and prune deleted remote references
git fetch -p

# Delete all local branches that no longer exist on remote
git fetch -p
git branch -vv | grep ': gone]' | awk '{print $1}' | xargs git branch -D