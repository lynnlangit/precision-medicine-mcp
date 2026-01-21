"""Base provider abstraction for LLM providers with MCP support.

Defines the interface that all providers must implement.
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Optional, Any
from dataclasses import dataclass


@dataclass
class ChatMessage:
    """Standardized chat message format."""
    role: str  # "user" or "assistant"
    content: str


@dataclass
class UsageInfo:
    """Token usage information."""
    input_tokens: int
    output_tokens: int
    total_tokens: int


@dataclass
class ChatResponse:
    """Standardized response format from providers."""
    content: str  # Main response text
    usage: Optional[UsageInfo] = None  # Token usage info
    raw_response: Any = None  # Provider-specific raw response


class LLMProvider(ABC):
    """Abstract base class for LLM providers with MCP support."""

    @abstractmethod
    def __init__(self, api_key: Optional[str] = None):
        """Initialize the provider.

        Args:
            api_key: API key for the provider (or None to use env var)
        """
        pass

    @abstractmethod
    def send_message(
        self,
        messages: List[ChatMessage],
        mcp_servers: List[Dict],
        model: str,
        max_tokens: int = 4096,
        temperature: float = 1.0,
        uploaded_files: Optional[Dict] = None
    ) -> ChatResponse:
        """Send message to the LLM with MCP servers enabled.

        Args:
            messages: Chat message history
            mcp_servers: MCP server configurations (format depends on provider)
            model: Model identifier
            max_tokens: Maximum response tokens
            temperature: Sampling temperature
            uploaded_files: Dict of uploaded files with paths and metadata

        Returns:
            ChatResponse with standardized format
        """
        pass

    @abstractmethod
    def get_provider_name(self) -> str:
        """Get the provider name for display.

        Returns:
            Provider name (e.g., "Claude", "Gemini")
        """
        pass

    @abstractmethod
    def get_model_display_name(self, model: str) -> str:
        """Get human-readable model name.

        Args:
            model: Model identifier

        Returns:
            Display name (e.g., "Claude Sonnet 4.5")
        """
        pass

    @abstractmethod
    def format_mcp_servers(self, server_configs: List[Dict]) -> List[Dict]:
        """Format MCP server configs for this provider's API.

        Different providers may expect different URL formats or structures.
        For example:
        - Anthropic expects: {type: "url", url: "...", name: "..."}
        - Gemini expects: {type: "mcp_server", name: "...", url: "..."}

        Args:
            server_configs: Standard MCP server configs from mcp_config.py

        Returns:
            Provider-specific formatted configs
        """
        pass

    def is_available(self) -> bool:
        """Check if provider is available (API key set, etc.).

        Returns:
            True if provider can be used
        """
        return True
