import argparse
import subprocess
import os


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="command"
    )

    args = parser.parse_args()

    token = subprocess.Popen(["gcloud", "auth", "print-access-token"],stdout=subprocess.PIPE).communicate()[0].decode("UTF-8")

    print(token)

    os.environ["GCS_OAUTH_TOKEN"] = token
    os.system(args.i)

