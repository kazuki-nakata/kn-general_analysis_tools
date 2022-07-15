# %%
import time
from hydra import compose, initialize_config_dir
import os
from datetime import datetime, timedelta

# %%


def hour_rounder(t):
    return (t.replace(second=0, microsecond=0, minute=0, hour=t.hour)
               +timedelta(hours=t.minute//30))
    
def deco_print_time(start_message: str, end_message: str):
    def _print_args(func):
        def new_function(*args, **kwargs):
            print(start_message, ": func name->", str(func))
            ftime = time.time()
            result = func(*args, **kwargs)
            ltime = time.time()
            print(end_message, ": time=", ltime - ftime)
            return result

        return new_function

    return _print_args


def import_config(config_path):
    conf_dir = os.path.dirname(os.path.abspath(config_path))
    conf_name = os.path.basename(config_path)
    with initialize_config_dir(config_dir=conf_dir):
        conf = compose(config_name=conf_name)
    return conf


def deco_import_config(config_path: str):
    def _perform_method(func):
        def new_function(*args):
            conf = import_config(config_path)
            print("used conf:", conf)
            result = func(*args, conf)
            return result

        return new_function

    return _perform_method


# %%
if __name__ == "__main__":

    @deco_print_time("処理開始", "処理終了")
    def test():
        print("test")
        return "test"

    output = test()

    # @deco_import_config(config_path=r"./dev/conf/test.yaml")
    # def test(a, b, c, conf):
    #     return a, b, c, conf

    output = test("a", "b", "c")
    print(output)
