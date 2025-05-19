from dotenv import load_dotenv
import sys
import os

load_dotenv(os.path.dirname(__file__) + os.sep + "../../.env")

GSHHG_PATH = os.environ.get("GSHHG_PATH")
PYRTTOV_PATH = os.environ.get("PYRTTOV_PATH")
sys.path.append(PYRTTOV_PATH)
