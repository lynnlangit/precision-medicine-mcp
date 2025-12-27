"""Retry utilities for external API calls with exponential backoff."""

import time
import logging
from typing import Callable, Any, Optional, Type, Tuple
from functools import wraps

logger = logging.getLogger(__name__)


def retry_with_backoff(
    max_retries: int = 3,
    base_delay: float = 1.0,
    max_delay: float = 30.0,
    exponential_base: float = 2.0,
    exceptions: Tuple[Type[Exception], ...] = (Exception,),
    on_retry: Optional[Callable[[Exception, int], None]] = None
):
    """
    Decorator for retrying functions with exponential backoff.

    Args:
        max_retries: Maximum number of retry attempts (default: 3)
        base_delay: Initial delay in seconds (default: 1.0)
        max_delay: Maximum delay between retries in seconds (default: 30.0)
        exponential_base: Base for exponential backoff (default: 2.0)
        exceptions: Tuple of exceptions to catch and retry (default: all Exception)
        on_retry: Optional callback function(exception, attempt_num) called before each retry

    Returns:
        Decorated function that retries on failure

    Example:
        ```python
        @retry_with_backoff(max_retries=3, exceptions=(requests.RequestException,))
        def fetch_tcga_data(cohort_id):
            response = requests.get(f"https://api.gdc.cancer.gov/cases?filters={cohort_id}")
            response.raise_for_status()
            return response.json()
        ```

    Retry schedule with default settings:
        - Attempt 1: Immediate
        - Attempt 2: Wait 1 second
        - Attempt 3: Wait 2 seconds
        - Attempt 4: Wait 4 seconds
        - (max 3 retries = 4 total attempts)
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def sync_wrapper(*args, **kwargs) -> Any:
            delay = base_delay
            last_exception = None

            for attempt in range(max_retries + 1):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    last_exception = e

                    if attempt == max_retries:
                        logger.error(
                            f"❌ {func.__name__} failed after {max_retries} retries: {e}"
                        )
                        raise

                    logger.warning(
                        f"⚠️  {func.__name__} attempt {attempt + 1}/{max_retries + 1} failed: {e}"
                    )

                    # Call on_retry callback if provided
                    if on_retry:
                        try:
                            on_retry(e, attempt + 1)
                        except Exception as callback_error:
                            logger.warning(f"on_retry callback failed: {callback_error}")

                    logger.info(f"⏳ Retrying in {delay:.1f}s...")
                    time.sleep(delay)

                    # Exponential backoff
                    delay = min(delay * exponential_base, max_delay)

            # Should never reach here, but just in case
            raise last_exception

        @wraps(func)
        async def async_wrapper(*args, **kwargs) -> Any:
            import asyncio
            delay = base_delay
            last_exception = None

            for attempt in range(max_retries + 1):
                try:
                    return await func(*args, **kwargs)
                except exceptions as e:
                    last_exception = e

                    if attempt == max_retries:
                        logger.error(
                            f"❌ {func.__name__} failed after {max_retries} retries: {e}"
                        )
                        raise

                    logger.warning(
                        f"⚠️  {func.__name__} attempt {attempt + 1}/{max_retries + 1} failed: {e}"
                    )

                    if on_retry:
                        try:
                            on_retry(e, attempt + 1)
                        except Exception as callback_error:
                            logger.warning(f"on_retry callback failed: {callback_error}")

                    logger.info(f"⏳ Retrying in {delay:.1f}s...")
                    await asyncio.sleep(delay)

                    delay = min(delay * exponential_base, max_delay)

            raise last_exception

        # Detect if function is async
        import asyncio
        if asyncio.iscoroutinefunction(func):
            return async_wrapper
        else:
            return sync_wrapper

    return decorator


def optional_api_call(
    fallback_value: Any = None,
    log_failure: bool = True,
    exceptions: Tuple[Type[Exception], ...] = (Exception,)
):
    """
    Decorator for optional API calls that should degrade gracefully.

    If the API call fails, returns fallback_value instead of raising exception.
    Use this for non-critical API calls where the analysis can continue without the data.

    Args:
        fallback_value: Value to return if API call fails (default: None)
        log_failure: Whether to log the failure as warning (default: True)
        exceptions: Tuple of exceptions to catch (default: all Exception)

    Returns:
        Decorated function that returns fallback_value on failure

    Example:
        ```python
        @optional_api_call(fallback_value=None)
        def fetch_tcga_cohort(cancer_type):
            # If TCGA API is down, return None instead of crashing
            return tcga_api.get_cohort(cancer_type)

        # Usage
        cohort = fetch_tcga_cohort("TCGA-OV")
        if cohort is None:
            logger.warning("TCGA API unavailable, skipping cohort comparison")
        else:
            # Process cohort data
            ...
        ```
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def sync_wrapper(*args, **kwargs) -> Any:
            try:
                return func(*args, **kwargs)
            except exceptions as e:
                if log_failure:
                    logger.warning(
                        f"⚠️  Optional API call {func.__name__} failed: {e}"
                    )
                    logger.info(
                        f"💡 Continuing with fallback value: {fallback_value}"
                    )
                return fallback_value

        @wraps(func)
        async def async_wrapper(*args, **kwargs) -> Any:
            try:
                return await func(*args, **kwargs)
            except exceptions as e:
                if log_failure:
                    logger.warning(
                        f"⚠️  Optional API call {func.__name__} failed: {e}"
                    )
                    logger.info(
                        f"💡 Continuing with fallback value: {fallback_value}"
                    )
                return fallback_value

        # Detect if function is async
        import asyncio
        if asyncio.iscoroutinefunction(func):
            return async_wrapper
        else:
            return sync_wrapper

    return decorator


class CircuitBreaker:
    """
    Circuit breaker pattern for API calls.

    Prevents repeated calls to a failing service by "opening the circuit"
    after a threshold of failures. After a timeout period, the circuit
    becomes "half-open" and allows one test call through.

    States:
        - CLOSED: Normal operation, requests pass through
        - OPEN: Service is failing, requests fail immediately
        - HALF_OPEN: Testing if service has recovered

    Example:
        ```python
        tcga_breaker = CircuitBreaker(
            failure_threshold=5,
            recovery_timeout=60.0,
            expected_exception=requests.RequestException
        )

        @tcga_breaker
        def fetch_tcga_data(cohort_id):
            response = requests.get(f"https://api.gdc.cancer.gov/...")
            response.raise_for_status()
            return response.json()
        ```
    """

    def __init__(
        self,
        failure_threshold: int = 5,
        recovery_timeout: float = 60.0,
        expected_exception: Type[Exception] = Exception,
        name: Optional[str] = None
    ):
        """
        Initialize circuit breaker.

        Args:
            failure_threshold: Number of failures before opening circuit
            recovery_timeout: Seconds to wait before attempting recovery
            expected_exception: Exception type that counts as failure
            name: Optional name for logging
        """
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.expected_exception = expected_exception
        self.name = name or "CircuitBreaker"

        self.failure_count = 0
        self.last_failure_time = None
        self.state = "CLOSED"  # CLOSED, OPEN, HALF_OPEN

    def __call__(self, func: Callable) -> Callable:
        """Decorator to apply circuit breaker to function."""
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            # Check if circuit is open
            if self.state == "OPEN":
                # Check if recovery timeout has passed
                if time.time() - self.last_failure_time >= self.recovery_timeout:
                    logger.info(
                        f"🔄 {self.name}: Attempting recovery (HALF_OPEN)"
                    )
                    self.state = "HALF_OPEN"
                else:
                    # Circuit is still open, fail fast
                    logger.warning(
                        f"⚠️  {self.name}: Circuit is OPEN, failing fast"
                    )
                    raise Exception(
                        f"Circuit breaker is OPEN for {func.__name__}. "
                        f"Service unavailable. Retry after "
                        f"{self.recovery_timeout - (time.time() - self.last_failure_time):.0f}s"
                    )

            # Attempt the call
            try:
                result = func(*args, **kwargs)

                # Success - reset circuit if it was half-open or had failures
                if self.state == "HALF_OPEN":
                    logger.info(f"✅ {self.name}: Service recovered (CLOSED)")
                    self.state = "CLOSED"
                    self.failure_count = 0

                return result

            except self.expected_exception as e:
                # Failure occurred
                self.failure_count += 1
                self.last_failure_time = time.time()

                logger.warning(
                    f"⚠️  {self.name}: Failure {self.failure_count}/{self.failure_threshold}"
                )

                # Check if we should open the circuit
                if self.failure_count >= self.failure_threshold:
                    self.state = "OPEN"
                    logger.error(
                        f"🔴 {self.name}: Circuit is now OPEN (too many failures)"
                    )

                # Re-raise the exception
                raise

        return wrapper

    def reset(self):
        """Manually reset the circuit breaker to CLOSED state."""
        logger.info(f"🔄 {self.name}: Manually reset to CLOSED")
        self.state = "CLOSED"
        self.failure_count = 0
        self.last_failure_time = None


# Example usage and testing
if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(level=logging.INFO)

    # Example 1: Retry with backoff
    @retry_with_backoff(max_retries=3, base_delay=0.5)
    def unreliable_function(success_rate: float = 0.3):
        import random
        if random.random() > success_rate:
            raise Exception("Random failure")
        return "Success!"

    # Example 2: Optional API call
    @optional_api_call(fallback_value={"status": "unavailable"})
    def fetch_external_data():
        raise Exception("API is down")

    # Example 3: Circuit breaker
    breaker = CircuitBreaker(failure_threshold=3, recovery_timeout=5.0)

    @breaker
    def failing_service():
        raise Exception("Service unavailable")

    # Test retry
    try:
        result = unreliable_function(success_rate=0.5)
        print(f"✅ Retry test result: {result}")
    except Exception as e:
        print(f"❌ Retry test failed: {e}")

    # Test optional API
    result = fetch_external_data()
    print(f"💡 Optional API result: {result}")

    # Test circuit breaker
    for i in range(10):
        try:
            failing_service()
        except Exception as e:
            print(f"Attempt {i+1}: {e}")
            time.sleep(0.5)
