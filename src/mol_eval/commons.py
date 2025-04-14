import functools
import json
import multiprocessing
import multiprocessing.pool
import os

from mol_eval.schemas import ConfigSchema
from mol_eval.logger import Logger

logging = Logger().get_logger()


def timeout(max_timeout):
    """Timeout decorator, parameter in seconds.

    Args:
        max_timeout (int): Maximum timeout in seconds.

    Returns:
        func_wrapper: The wrapped function.
    """

    def timeout_decorator(item):
        """
        Wrap the original function.

        Args:
            item: The function to wrap.

        Returns:
            func_wrapper: The wrapped function.
        """

        @functools.wraps(item)
        def func_wrapper(*args, **kwargs):
            """Closure for function."""
            pool = multiprocessing.pool.ThreadPool(processes=1)
            async_result = pool.apply_async(item, args, kwargs)
            try:
                # Try to get the result within the max_timeout
                return async_result.get(timeout=max_timeout)
            except multiprocessing.TimeoutError:
                # Catch the TimeoutError from multiprocessing
                logging.warning(
                    f"Function {item.__name__} exceeded the timeout of {max_timeout} seconds."
                )
                return []  # Return empty list or any default value to indicate timeout
            except Exception as e:
                # Catch any other exceptions
                logging.error(f"An error occurred: {e}")
                return []  # Gracefully handle other exceptions
            finally:
                pool.close()

        return func_wrapper

    return timeout_decorator


def load_config_file(file_path: str, read_op_type: str = "r") -> str:
    """Load a file and return its content.

    Args:
        file_path (str): Path to the file.
        read_op_type (str): Type of read operation (e.g., 'r' for read, 'rb' for binary read).
    Returns:
        str: Content of the file.
    """
    if file_path.split(".")[-1] != "json":
        raise ValueError("File must be a JSON file.")
    if not file_path:
        raise ValueError("File path cannot be None or empty.")
    if not read_op_type:
        raise Warning("Read operation type is not specified. Defaulting to 'r'.")
    if read_op_type not in ["r", "rb"]:
        raise ValueError(
            "Invalid read operation type. Use 'r' for text or 'rb' for binary."
        )
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    with open(file_path, read_op_type) as file:
        config_data = json.load(file)

    return ConfigSchema(**config_data)
