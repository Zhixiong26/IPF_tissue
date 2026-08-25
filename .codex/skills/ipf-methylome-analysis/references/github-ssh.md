# GitHub and SSH setup

Use this reference when a repository has not yet been provided, the local project is not connected, or SSH authentication fails.

## Information to request

- GitHub SSH clone URL, normally `git@github.com:OWNER/REPO.git`.
- Intended default branch, normally `main` unless the user specifies another branch.
- Whether the repository already exists and whether it contains commits.
- Whether the user wants only local setup or also authorizes an initial commit/push.

Never ask for or display a private SSH key, GitHub password, personal access token, recovery code, or session cookie. A public key ending in `.pub`, its fingerprint, and non-secret repository metadata are safe to provide.

## Safe verification sequence

1. Inspect `git status`, the current branch, and `git remote -v` without changing them.
2. Check whether an SSH public key exists without reading any private-key contents.
3. Test GitHub authentication with `ssh -T git@github.com` only when network access and the user's environment permit it. This test authenticates but does not push repository data.
4. If no key exists, propose an Ed25519 key tied to the user's account email. Key creation changes the user's SSH configuration, so obtain confirmation immediately before doing it.
5. Ask the user to add the resulting `.pub` key to GitHub, then repeat authentication testing.
6. Add or change `origin` only after confirming the exact URL. Fetch before reconciling a non-empty remote with a non-empty local repository.

## Repository hygiene

Before the initial commit, inspect or create `.gitignore` rules appropriate to the project. At minimum keep these classes out of Git:

- `Data/Raw_fastq/` and other linked or copied raw sequencing data.
- Large `Results/` products, alignment files, coverage matrices, model checkpoints, and generated figures unless the user intentionally selects small deliverables.
- `Scripts/*/logs/`, legacy root-level `logs/` and `log/`, scheduler stdout/stderr, caches, environments, credentials, and private keys.

Keep scripts, small configuration files, sample-sheet schemas, maintained supplementary metadata, documentation, and the project skill versioned when they contain no sensitive data.
