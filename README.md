### Course material for Single Cell Analysis in Python 2024, EBI

[Reference](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#create-and-use-virtual-environments)

Unix/macOS

1. Install Python [v3.9.6](https://www.python.org/downloads/release/python-396/) and check installation `python3 --version` should return `Python 3.9.6`
2. Clone this repo and enter this repo `cd <path_to_this_repo>`
3. Create a virtual environment by ways of choice, for instance, `python3 -m venv .venv`
4. Activate the virtial environment by `source .venv/bin/activate`
5. Confirm the correct python is being used `which python` should return `<path_to_this_repo>/.venv/bin/python` and `python --version` should return `Python 3.9.6`
6. Prepare pip `python3 -m pip install --upgrade pip`
7. Install the required packages by `pip install -r requirements.txt`

credit:

- materials: Hugo Tavares, Universiy of Cambridge, UK
- conversions: Diana Yu, EMBL-EBI, UK
- envs: Yuyao Song, EMBL-EBI, UK
