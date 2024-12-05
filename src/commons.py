import functools
import multiprocessing
import multiprocessing.pool


def timeout(max_timeout):
    """
    Timeout decorator, parameter in seconds.

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
                print(f"Function {item.__name__} exceeded the timeout of {max_timeout} seconds.")
                return []  # Return empty list or any default value to indicate timeout
            except Exception as e:
                # Catch any other exceptions
                print(f"An error occurred: {e}")
                return []  # Gracefully handle other exceptions
            finally:
                pool.close()

        return func_wrapper

    return timeout_decorator
