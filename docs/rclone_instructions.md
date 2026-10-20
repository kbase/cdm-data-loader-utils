# `rclone` Instructions

`rclone` is a tool for syncing files across various cloud and local storage systems. Detailed descriptions
of installation, usage, and supported protocols can be found in their documentation: https://rclone.org/

Included here are simplified instructions for downloading and configuring `rclone` in the lakehouse
environment, and using it to sync a few records from PDB (an `http` remote) as an example.

## Download and Install

Start a lakehouse session and open a terminal window. Then, download `rclone` and put it somewhere that's
in your `PATH`. I just created a `bin/` folder in my home directory.
```sh
$ curl -O https://downloads.rclone.org/rclone-current-linux-amd64.zip
$ bsdtar -xf rclone-current-linux-amd64.zip
$ mkdir ~/bin
$ cp rclone-v1.75.0-linux-amd64/rclone ~/bin/
$ cd ~/bin
$ chmod 755 rclone
$ export PATH=$(pwd):$PATH
```

You'll probably want to modify your `.custom_profile` to include `~/bin` on your `PATH` in every
new `bash` session.

## Configure Remotes

Configuring `rclone` remotes is an interactive process started by running:

```sh
$ rclone config
```

You'll see any already configured remotes, and options for adding/removing/modifying remotes:

```sh
$ rclone config
Current remotes:

Name                 Type
====                 ====

e) Edit existing remote
n) New remote
d) Delete remote
r) Rename remote
c) Copy remote
s) Set configuration password
q) Quit config
e/n/d/r/c/s/q> 
```

Detailed options for each remote type can be found here: https://rclone.org/overview/

We'll create two remotes: one for the lakehouse S3 store (type: `AWS S3: Minio`) and
one for PDB (type: `HTTP`)

### Lakehouse Remote

This remote can be used as the destination for all downloads. From the `rclone config`
main page, choose `n` (New remote), then enter the following when prompted (replacing
`my-access-key-id` and `my-secret-access-key` with your actual lakehouse credentials). There are many options listed for several steps. I'm only showing a few options
including the one you should choose from the list.

```
e/n/d/r/c/s/q> n
name> lakehouse

Option Storage.
Type of storage to configure.
Choose a number from below, or type in your own value.
 4 / Amazon S3 Compliant Storage Providers including AWS, Alibaba, ArvanCloud, BizflyCloud, Ceph, ChinaMobile, Cloudflare, Cubbit, DigitalOcean, Dreamhost, Exaba, Fastly, FileLu, FlashBlade, GCS, HCP, Hetzner, HuaweiOBS, IBMCOS, IDrive, ImpossibleCloud, Intercolo, IONOS, Leviia, Liara, Linode, LyveCloud, Magalu, Mega, Minio, Netease, Outscale, OVHcloud, Petabox, Qiniu, Rabata, RackCorp, Rclone, Scaleway, Scality, SeaweedFS, Selectel, Servercore, SpectraLogic, Storj, Synology, TencentCOS, US3, Wasabi, Zadara, Zata, ZeroServices, Other - Tier 1
Storage> 4

Option provider.
Choose your S3 provider.
Choose a number from below, or type in your own value.
Press Enter to leave empty.
 5 / Ceph Object Storage
   \ (Ceph)
 30 / Minio Object Storage
   \ (Minio)
provider> 30

Option env_auth.
Get AWS credentials from runtime (environment variables or EC2/ECS meta data if no env vars).
Only applies if access_key_id and secret_access_key is blank.
Choose a number from below, or type in your own boolean value (true or false).
Press Enter for the default (false).
 1 / Enter AWS credentials in the next step.
   \ (false)
 2 / Get AWS credentials from the environment (env vars or IAM).
   \ (true)
env_auth> 1

Option access_key_id.
AWS Access Key ID.
Leave blank for anonymous access or runtime credentials.
Enter a value. Press Enter to leave empty.
access_key_id> my-access-key-id

Option secret_access_key.
AWS Secret Access Key (password).
Leave blank for anonymous access or runtime credentials.
Enter a value. Press Enter to leave empty.
secret_access_key> my-secret-access-key

Option region.
Region to connect to.
Leave blank if you are using an S3 clone and you don't have a region.
Choose a number from below, or type in your own value.
Press Enter to leave empty.
   / Use this if unsure.
 1 | Will use v4 signatures and an empty region.
   \ ()
   / Use this only if v4 signatures don't work.
 2 | E.g. pre Jewel/v10 CEPH.
   \ (other-v2-signature)
region> 

Option endpoint.
Endpoint for S3 API.
Required when using an S3 clone.
Enter a value. Press Enter to leave empty.
endpoint> https://minio.berdl.kbase.us

Option location_constraint.
Location constraint - must be set to match the Region.
Leave blank if not sure. Used when creating buckets only.
Enter a value. Press Enter to leave empty.
location_constraint> 

Option acl.
Canned ACL used when creating buckets and storing or copying objects.
This ACL is used for creating objects and if bucket_acl isn't set, for creating buckets too.
For more info visit https://docs.aws.amazon.com/AmazonS3/latest/dev/acl-overview.html#canned-acl
Note that this ACL is applied when server-side copying objects as S3
doesn't copy the ACL from the source but rather writes a fresh one.
If the acl is an empty string then no X-Amz-Acl: header is added and
the default (private) will be used.
Choose a number from below, or type in your own value.
Press Enter to leave empty.
   / Owner gets FULL_CONTROL.
 1 | No one else has access rights (default).
   \ (private)
   / Owner gets FULL_CONTROL.
 2 | The AllUsers group gets READ access.
   \ (public-read)
acl> 2

Option server_side_encryption.
The server-side encryption algorithm used when storing this object in S3.
Choose a number from below, or type in your own value.
Press Enter to leave empty.
 1 / None
   \ ()
 2 / AES256
   \ (AES256)
 3 / aws:kms
   \ (aws:kms)
server_side_encryption> 

Option sse_kms_key_id.
If using KMS ID you must provide the ARN of Key.
Choose a number from below, or type in your own value.
Press Enter to leave empty.
 1 / None
   \ ()
 2 / arn:aws:kms:*
   \ (arn:aws:kms:us-east-1:*)
sse_kms_key_id> 

Option bucket_object_lock_enabled.
Enable Object Lock when creating new buckets.
Enter a boolean value (true or false). Press Enter for the default (false).
bucket_object_lock_enabled> 

Edit advanced config?
y) Yes
n) No (default)
y/n> n

Configuration complete.
Options:
- type: s3
- provider: Minio
- access_key_id: my-access-key-id
- secret_access_key: my-secret-access-key
- endpoint: https://minio.berdl.kbase.us
- acl: public-read
Keep this "lakehorse" remote?
y) Yes this is OK (default)
e) Edit this remote
d) Delete this remote
y/e/d> y
```

You could select `2` for `env_auth` if your lakehouse S3 credentials are store in the standard AWS
environment variables.

You also might want another option for `acl` depending on what permissions you want for your
transferred files.

### PDB Reomte

Next, let's configure a remote for PDB (and https/ftp endpoint). From the main page of `rclone config`
choose `n` (New remote) and enter the following information:

```
e/n/d/r/c/s/q> n
name> pdb

Option Storage.
Type of storage to configure.
Choose a number from below, or type in your own value.
26 / HTTP - Tier 3
   \ (http)
Storage> 26

Option url.
URL of HTTP host to connect to.
E.g. "https://example.com", or "https://user:pass@example.com" to use a username and password.
Enter a value.
url> https://ftp.ebi.ac.uk

Option no_escape.
Do not escape URL metacharacters in path names.
Enter a boolean value (true or false). Press Enter for the default (false).
no_escape> 

Edit advanced config?
y) Yes
n) No (default)
y/n> 

Configuration complete.
Options:
- type: http
- url: https://ftp.ebi.ac.uk
Keep this "pdb" remote?
y) Yes this is OK (default)
e) Edit this remote
d) Delete this remote
y/e/d> y
```

Super simple. Now select `q` to exit `rclone config`.

## Transfer

Now let's transfer some files. I'm going to transfer a small set of PDB records to my personal S3 "folder".

First, run the sync with `--dry-run` to see which files will be transferred and how big they are:

```sh
$ rclone sync --dry-run pdb-https:pub/databases/wwpdb/pdb/data/entries/06/ lakehouse:cdm-lake/users-general-warehouse/mattdawson/pdb/raw_data/06/
...
2026/08/19 19:09:52 NOTICE: pdb_0000306d/validation_reports/pdb_0000306d_validation.cif.gz: Skipped copy as --dry-run is set (size 2.406Ki)
2026/08/19 19:09:52 NOTICE: pdb_0000306d/validation_reports/pdb_0000306d_validation.pdf.gz: Skipped copy as --dry-run is set (size 397.341Ki)
2026/08/19 19:09:52 NOTICE: pdb_0000306d/validation_reports/pdb_0000306d_multipercentile_validation.svg.gz: Skipped copy as --dry-run is set (size 651)
2026/08/19 19:09:52 NOTICE: pdb_0000306d/validation_reports/pdb_0000306d_validation.xml.gz: Skipped copy as --dry-run is set (size 1.550Ki)
2026/08/19 19:09:52 ERROR : Attempt 2/3 succeeded
2026/08/19 19:09:52 NOTICE: 
Transferred:       10.919 MiB / 10.919 MiB, 100%, 0 B/s, ETA -
Checks:                 0 / 0, -, Listed 149
Transferred:          101 / 101, 100%
Elapsed time:         2.7s
$
```

Then, to do the actual transfer, just remove the `--dry-run` and do it again:

```sh
$ rclone sync pdb-https:pub/databases/wwpdb/pdb/data/entries/06/ lakehouse:cdm-lake/users-general-warehouse/mattdawson/pdb/raw_data/06/
```

You can check that the files were transferred with an `ls` of the lakehouse folder:

```sh
$ rclone ls lakehouse:cdm-lake/users-general-warehouse/mattdawson/pdb/raw_data/06/
     9641 pdb_0000106d/assemblies/pdb_0000106d-assembly1.cif.gz
    53000 pdb_0000106d/structures/pdb_0000106d-extatom.xml.gz
    10202 pdb_0000106d/structures/pdb_0000106d-noatom.xml.gz
    66717 pdb_0000106d/structures/pdb_0000106d.cif.gz
    41342 pdb_0000106d/structures/pdb_0000106d.pdb.gz
    88796 pdb_0000106d/structures/pdb_0000106d.xml.gz
...
$
```