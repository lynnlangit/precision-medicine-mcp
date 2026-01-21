"""LLM Provider abstraction for MCP chat.

Provides a unified interface for different LLM providers (Claude, Gemini).
"""

from typing import Optional
import os

from .base import LLMProvider, ChatMessage, ChatResponse, UsageInfo
from .anthropic_provider import AnthropicProvider

# Import Gemini provider conditionally
try:
    from .gemini_provider import GeminiProvider
    GEMINI_AVAILABLE = True
except ImportError:
    GeminiProvider = None
    GEMINI_AVAILABLE = False


# Export public API
__all__ = [
    "LLMProvider",
    "ChatMessage",
    "ChatResponse",
    "UsageInfo",
    "AnthropicProvider",
    "GeminiProvider",
    "get_provider",
    "get_available_providers",
    "GEMINI_AVAILABLE"
]


def get_provider(provider_name: str = "claude", api_key: Optional[str] = None) -> LLMProvider:
    """Factory function to get a provider instance.

    Args:
        provider_name: Provider name ("claude" or "gemini")
        api_key: Optional API key (uses env var if not provided)

    Returns:
        Provider instance

    Raises:
        ValueError: If provider not found or not available
    """
    provider_name = provider_name.lower()

    if provider_name == "claude":
        return AnthropicProvider(api_key=api_key)

    elif provider_name == "gemini":
        if not GEMINI_AVAILABLE or GeminiProvider is None:
            raise ValueError(
                "Gemini provider not available. "
                "Install with: pip install google-genai>=1.55.0"
            )
        return GeminiProvider(api_key=api_key)

    else:
        raise ValueError(f"Unknown provider: {provider_name}")


def get_available_providers() -> dict:
    """Get list of available providers with their status.

    Returns:
        Dict mapping provider names to availability info
    """
    providers = {}

    # Check Claude
    claude_key = os.getenv("ANTHROPIC_API_KEY")
    providers["claude"] = {
        "name": "Claude",
        "available": bool(claude_key),
        "models": [
            {"id": "claude-sonnet-4-20250514", "name": "Claude Sonnet 4.5"},
            {"id": "claude-opus-4-5", "name": "Claude Opus 4.5"},
            {"id": "claude-haiku-4", "name": "Claude Haiku 4"}
        ],
        "default_model": "claude-sonnet-4-20250514"
    }

    # Check Gemini
    gemini_key = os.getenv("GEMINI_API_KEY")
    providers["gemini"] = {
        "name": "Gemini",
        "available": GEMINI_AVAILABLE and bool(gemini_key),
        "models": [
            {"id": "gemini-3-flash-preview", "name": "Gemini 3.0 Flash"}
        ],
        "default_model": "gemini-3-flash-preview",
        "note": "Gemini API key required" if GEMINI_AVAILABLE else "Install google-genai package"
    }

    return providers
