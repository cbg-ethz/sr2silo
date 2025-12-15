# Multi-Virus Deployments

This directory contains the configuration and scripts for deploying the sr2silo workflow across multiple viruses (e.g., SARS-CoV-2, RSV-A) on the ETH Euler cluster.

## Structure

```
deployments/
├── secrets.py            # Credential loader (Python) - decrypts secrets/credentials.enc
├── secrets/              # Shared Loculus credentials (GPG encrypted)
│   └── credentials.enc   # Encrypted file containing USERNAME and PASSWORD
├── submit-daily.sbatch   # Shared SLURM script for daily submission and auto-resubmission
├── _templates/           # Template directory for adding new viruses
├── covid/                # SARS-CoV-2 deployment
│   └── config.yaml       # Virus-specific configuration
└── rsva/                 # RSV-A deployment
    └── config.yaml       # Virus-specific configuration
```

## 1. Setup Credentials (One-time)

We use GPG symmetric encryption to store the shared Loculus credentials securely. You need to know the passphrase to encrypt/decrypt them.

**To create or update the encrypted credentials:**

```bash
# You will be prompted for the passphrase twice
echo -e "USERNAME=vpipe\nPASSWORD=your_actual_password" | gpg --symmetric --cipher-algo AES256 -o deployments/secrets/credentials.enc
```

## 2. Run Deployment

The deployment uses a single SLURM script (`submit-daily.sbatch`) that takes the virus name as an argument. It runs the workflow and **automatically resubmits itself** for the next day at 23:00.

### Start a Daily Workflow

To start the daily cycle for a specific virus (e.g., `rsva`):

```bash
# Submit immediately (starts now)
sbatch --job-name=sr2silo-rsva --export=VIRUS=rsva,GPG_PASSPHRASE="your_passphrase" --begin=now deployments/submit-daily.sbatch

# Or schedule for tonight at 23:00
sbatch --job-name=sr2silo-rsva --export=VIRUS=rsva,GPG_PASSPHRASE="your_passphrase" --begin=23:00 deployments/submit-daily.sbatch
```

**Note:** The `GPG_PASSPHRASE` is required to decrypt the credentials at runtime. It is passed securely to the job environment.

### Monitor Jobs

The jobs will appear in the queue with virus-specific names (e.g., `sr2silo-rsva`, `sr2silo-covid`).

```bash
squeue -u $USER -o "%.18i %.9P %.20j %.8u %.2t %.10M %.6D %R"
```

### Stop the Cycle

To stop the daily resubmission, simply cancel the pending job:

```bash
scancel <JOB_ID>
```

## 3. Add a New Virus

1.  Copy the template directory:
    ```bash
    cp -r deployments/_templates deployments/my-new-virus
    ```
2.  Edit `deployments/my-new-virus/config.yaml`:
    *   Update `VIRUS_NAME` and `VIRUS_DISPLAY_NAME`.
    *   Set the correct `ORGANISM` (must match sr2silo supported organisms).
    *   Update `LOCATIONS`, `BASE_SAMPLE_DIR`, and `TIMELINE_FILE` paths.
    *   Adjust `GROUP_ID` if needed.
3.  Test the configuration:
    ```bash
    cd workflow
    snakemake --configfile ../deployments/my-new-virus/config.yaml -n
    ```
4.  Start the deployment as described in step 2.

