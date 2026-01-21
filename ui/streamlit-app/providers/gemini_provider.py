"""Google Gemini provider implementation.

Uses Google's Interactions API with remote MCP server support.
"""

import os
import asyncio
from typing import List, Dict, Optional, Any

try:
    from google import genai
    from google.genai import types
    GEMINI_AVAILABLE = True
except ImportError:
    GEMINI_AVAILABLE = False
    genai = None
    types = None

from .base import LLMProvider, ChatMessage, ChatResponse, UsageInfo


class GeminiProvider(LLMProvider):
    """Google Gemini provider with MCP support via Interactions API."""

    def __init__(self, api_key: Optional[str] = None):
        """Initialize Gemini provider.

        Args:
            api_key: Google AI API key (defaults to GEMINI_API_KEY env var)
        """
        if not GEMINI_AVAILABLE:
            raise ImportError(
                "google-genai package not installed. "
                "Install with: pip install google-genai>=1.55.0"
            )

        self.api_key = api_key or os.getenv("GEMINI_API_KEY")
        if not self.api_key:
            raise ValueError("GEMINI_API_KEY not found in environment")

        self.client = genai.Client(api_key=self.api_key)

    def send_message(
        self,
        messages: List[ChatMessage],
        mcp_servers: List[Dict],
        model: str = "gemini-3-flash-preview",
        max_tokens: int = 4096,
        temperature: float = 1.0,
        uploaded_files: Optional[Dict] = None
    ) -> ChatResponse:
        """Send message to Gemini with MCP servers.

        Args:
            messages: Chat message history
            mcp_servers: MCP server configurations
            model: Gemini model to use
            max_tokens: Maximum response tokens
            temperature: Sampling temperature
            uploaded_files: Dict of uploaded files

        Returns:
            ChatResponse with standardized format
        """
        # Format MCP servers for Gemini
        formatted_servers = self.format_mcp_servers(mcp_servers)

        # Convert messages to Gemini format
        input_messages = self._convert_messages_to_gemini_format(messages)

        # Build system prompt
        system_prompt = self._build_system_prompt(formatted_servers, uploaded_files)

        # Add system prompt as first message
        if system_prompt:
            input_messages.insert(0, {
                "role": "user",
                "content": [{"type": "text", "text": system_prompt}]
            })
            # Add a system-style response
            input_messages.insert(1, {
                "role": "model",
                "content": [{"type": "text", "text": "Understood. I will use the available MCP tools to perform actual analyses."}]
            })

        try:
            # Create config with thinking level
            config = types.GenerateContentConfig(
                temperature=temperature,
                max_output_tokens=max_tokens,
                thinking_level="medium"  # Gemini 3 specific parameter
            )

            # Use asyncio to run the async interaction
            interaction = asyncio.run(
                self._create_interaction_async(
                    model=model,
                    input_messages=input_messages,
                    tools=formatted_servers,
                    config=config
                )
            )

            # Extract content from interaction outputs
            content = self._format_response(interaction)

            # Extract usage info (Gemini may not provide this in the same way)
            usage = self._get_usage_info(interaction)

            return ChatResponse(
                content=content,
                usage=usage,
                raw_response=interaction
            )

        except Exception as e:
            raise Exception(f"Gemini API error: {str(e)}")

    async def _create_interaction_async(
        self,
        model: str,
        input_messages: List[Dict],
        tools: List[Dict],
        config: Any
    ):
        """Create interaction with Gemini API (async).

        Args:
            model: Model identifier
            input_messages: Formatted messages
            tools: MCP server tools
            config: Generation config

        Returns:
            Interaction response
        """
        return self.client.interactions.create(
            model=model,
            input=input_messages,
            tools=tools,
            config=config
        )

    def get_provider_name(self) -> str:
        """Get provider name."""
        return "Gemini"

    def get_model_display_name(self, model: str) -> str:
        """Get human-readable model name."""
        model_names = {
            "gemini-3-flash-preview": "Gemini 3.0 Flash",
            "gemini-3.0-flash-preview": "Gemini 3.0 Flash",
        }
        return model_names.get(model, model)

    def format_mcp_servers(self, server_configs: List[Dict]) -> List[Dict]:
        """Format MCP servers for Gemini Interactions API.

        Gemini expects:
        {
            "type": "mcp_server",
            "name": "server_name",
            "url": "https://server/sse"
        }

        Args:
            server_configs: Standard MCP server configs

        Returns:
            Gemini-formatted server configs
        """
        formatted = []
        for server in server_configs:
            url = server["url"]
            # Ensure URL has /sse suffix
            if not url.endswith("/sse"):
                url = url + "/sse"

            formatted.append({
                "type": "mcp_server",
                "name": server["name"],
                "url": url
            })

        return formatted

    def is_available(self) -> bool:
        """Check if Gemini provider is available."""
        return GEMINI_AVAILABLE and bool(self.api_key)

    def _convert_messages_to_gemini_format(self, messages: List[ChatMessage]) -> List[Dict]:
        """Convert ChatMessage objects to Gemini format.

        Gemini expects:
        [
            {"role": "user", "content": [{"type": "text", "text": "..."}]},
            {"role": "model", "content": [{"type": "text", "text": "..."}]}
        ]

        Args:
            messages: List of ChatMessage objects

        Returns:
            Gemini-formatted messages
        """
        gemini_messages = []
        for msg in messages:
            # Convert "assistant" role to "model" for Gemini
            role = "model" if msg.role == "assistant" else msg.role

            gemini_messages.append({
                "role": role,
                "content": [{"type": "text", "text": msg.content}]
            })

        return gemini_messages

    def _build_system_prompt(
        self,
        mcp_servers: List[Dict],
        uploaded_files: Optional[Dict] = None
    ) -> str:
        """Build system prompt for Gemini."""
        server_descriptions = [
            f"- {server['name']}: Available for bioinformatics analysis"
            for server in mcp_servers
        ]

        system_prompt = f"""You are a bioinformatics assistant with access to specialized analysis tools via MCP servers.

IMPORTANT: When users ask you to perform analysis, you MUST use the appropriate MCP tool rather than providing theoretical responses.

Available MCP servers:
{chr(10).join(server_descriptions)}

Guidelines:
- When asked to perform pathway enrichment, spatial analysis, or other bioinformatics tasks, USE THE TOOLS
- Call the appropriate MCP tool to get real analysis results
- After receiving tool results, interpret and explain them to the user
- Provide actionable insights based on the actual analysis output
- If you don't have access to a required tool, let the user know

Do not simulate or describe what an analysis would show - actually perform it using the available tools."""

        # Add file information if present (similar to Anthropic)
        if uploaded_files and len(uploaded_files) > 0:
            file_count = len(uploaded_files)
            file_list = []
            for filename, file_info in uploaded_files.items():
                original_name = file_info.get('original_name', filename)
                file_list.append(f"- {original_name}")

            system_prompt += f"""

UPLOADED FILES:
The user has uploaded {file_count} file(s): {', '.join(file_list)}
When analyzing these files, use the appropriate MCP tools with the file paths provided."""

        return system_prompt

    def _format_response(self, interaction) -> str:
        """Format Gemini response for display.

        Args:
            interaction: Gemini interaction response

        Returns:
            Formatted text response
        """
        if not hasattr(interaction, 'outputs') or not interaction.outputs:
            return "No response from Gemini"

        parts = []
        for output in interaction.outputs:
            # Extract text from output
            if hasattr(output, 'text'):
                parts.append(output.text)
            elif hasattr(output, 'content'):
                # Handle structured content
                for content_part in output.content:
                    if hasattr(content_part, 'text'):
                        parts.append(content_part.text)

        return "\n\n".join(parts) if parts else "Empty response"

    def _get_usage_info(self, interaction) -> Optional[UsageInfo]:
        """Extract usage information from Gemini response.

        Note: Gemini's usage reporting may differ from Anthropic.
        This is a best-effort extraction.

        Args:
            interaction: Gemini interaction response

        Returns:
            UsageInfo if available, None otherwise
        """
        # Gemini may provide usage info differently
        # Check for usage metadata in the interaction
        if hasattr(interaction, 'metadata') and hasattr(interaction.metadata, 'usage'):
            usage = interaction.metadata.usage
            if hasattr(usage, 'prompt_token_count') and hasattr(usage, 'candidates_token_count'):
                return UsageInfo(
                    input_tokens=usage.prompt_token_count,
                    output_tokens=usage.candidates_token_count,
                    total_tokens=usage.prompt_token_count + usage.candidates_token_count
                )

        # If usage info not available, return None
        return None
