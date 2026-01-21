"""Anthropic (Claude) provider implementation.

Wraps existing Claude API logic to maintain backward compatibility.
"""

import os
from typing import List, Dict, Optional
import anthropic

from .base import LLMProvider, ChatMessage, ChatResponse, UsageInfo


class AnthropicProvider(LLMProvider):
    """Anthropic Claude provider with MCP support."""

    def __init__(self, api_key: Optional[str] = None):
        """Initialize Anthropic provider.

        Args:
            api_key: Anthropic API key (defaults to ANTHROPIC_API_KEY env var)
        """
        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError("ANTHROPIC_API_KEY not found in environment")

        self.client = anthropic.Anthropic(api_key=self.api_key)

    def send_message(
        self,
        messages: List[ChatMessage],
        mcp_servers: List[Dict],
        model: str = "claude-sonnet-4-20250514",
        max_tokens: int = 4096,
        temperature: float = 1.0,
        uploaded_files: Optional[Dict] = None
    ) -> ChatResponse:
        """Send message to Claude with MCP servers.

        Args:
            messages: Chat message history
            mcp_servers: MCP server configurations
            model: Claude model to use
            max_tokens: Maximum response tokens
            temperature: Sampling temperature
            uploaded_files: Dict of uploaded files

        Returns:
            ChatResponse with standardized format
        """
        # Convert ChatMessage objects to dict format for Anthropic API
        api_messages = [
            {"role": msg.role, "content": msg.content}
            for msg in messages
        ]

        # Format MCP servers for Anthropic (already in correct format)
        formatted_servers = self.format_mcp_servers(mcp_servers)

        # Get tools config
        tools = [
            {"type": "mcp_toolset", "mcp_server_name": server["name"]}
            for server in formatted_servers
        ]

        # Build system prompt
        system_prompt = self._build_system_prompt(formatted_servers, uploaded_files)

        try:
            response = self.client.beta.messages.create(
                model=model,
                max_tokens=max_tokens,
                temperature=temperature,
                messages=api_messages,
                mcp_servers=formatted_servers,
                tools=tools,
                system=system_prompt,
                betas=["mcp-client-2025-11-20"]
            )

            # Extract text from response
            content = self._format_response(response)

            # Extract usage info
            usage = self._get_usage_info(response)

            return ChatResponse(
                content=content,
                usage=usage,
                raw_response=response
            )

        except anthropic.APIError as e:
            raise Exception(f"Claude API error: {str(e)}")
        except Exception as e:
            raise Exception(f"Unexpected error: {str(e)}")

    def get_provider_name(self) -> str:
        """Get provider name."""
        return "Claude"

    def get_model_display_name(self, model: str) -> str:
        """Get human-readable model name."""
        model_names = {
            "claude-sonnet-4-20250514": "Claude Sonnet 4.5",
            "claude-sonnet-4-5": "Claude Sonnet 4.5",
            "claude-opus-4-5": "Claude Opus 4.5",
            "claude-haiku-4": "Claude Haiku 4"
        }
        return model_names.get(model, model)

    def format_mcp_servers(self, server_configs: List[Dict]) -> List[Dict]:
        """Format MCP servers for Anthropic API.

        Anthropic expects:
        {
            "type": "url",
            "url": "https://server/sse",
            "name": "server_name"
        }

        Args:
            server_configs: Standard MCP server configs

        Returns:
            Anthropic-formatted server configs
        """
        formatted = []
        for server in server_configs:
            url = server["url"]
            # Ensure URL has /sse suffix for Anthropic
            if not url.endswith("/sse"):
                url = url + "/sse"

            formatted.append({
                "type": "url",
                "url": url,
                "name": server["name"]
            })

        return formatted

    def is_available(self) -> bool:
        """Check if Claude provider is available."""
        return bool(self.api_key)

    def _build_system_prompt(
        self,
        mcp_servers: List[Dict],
        uploaded_files: Optional[Dict] = None
    ) -> str:
        """Build system prompt with MCP server and file information."""
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

        # Add uploaded file information if present
        if uploaded_files and len(uploaded_files) > 0:
            file_descriptions = []
            for filename, file_info in uploaded_files.items():
                metadata = file_info['metadata']
                path = file_info['path']
                original_name = file_info.get('original_name', filename)
                is_gcs = file_info.get('source') == 'gcs'

                if is_gcs:
                    file_desc = f"""- **{original_name}** (GCS)
  - GCS URI: {path}
  - Bucket: {metadata.get('bucket', 'unknown')}
  - Path: {metadata.get('blob_path', 'unknown')}
  - Type: {metadata['extension']}"""
                else:
                    file_desc = f"""- **{original_name}**
  - File path: {path}
  - Type: {metadata['extension']}
  - Size: {metadata['size_mb']:.2f} MB"""

                if not metadata.get('is_binary', False):
                    file_desc += f"\n  - Lines: {metadata.get('line_count', 'N/A')}"

                    # Include content preview for small files
                    content = metadata.get('_gcs_content')
                    if not content and not is_gcs and metadata.get('size_bytes', 0) < 50000:
                        try:
                            with open(path, 'r') as f:
                                content = f.read()
                        except Exception:
                            pass

                    if content:
                        if len(content) > 10000:
                            content = content[:10000] + "\n\n[... truncated ...]"
                        file_desc += f"\n  - Content preview:\n```\n{content}\n```"

                file_descriptions.append(file_desc)

            system_prompt += f"""

UPLOADED FILES AVAILABLE:
The user has uploaded {len(uploaded_files)} file(s) for analysis.

{chr(10).join(file_descriptions)}

FILE ACCESS GUIDELINES:
1. **Small text files with content preview**: Analyze the content directly without calling MCP tools
2. **GCS files (gs://...)**: MCP servers on Cloud Run can access these directly - pass the GCS URI to MCP tool calls
3. **Local file paths**: MCP servers may not access these - prefer analyzing inline content

When calling MCP tools with GCS files, use the GCS URI (gs://bucket/path) as the file path parameter."""

        return system_prompt

    def _format_response(self, response: anthropic.types.Message) -> str:
        """Format Claude response for display."""
        if not response.content:
            return "No response from Claude"

        parts = []
        for block in response.content:
            if hasattr(block, 'text'):
                parts.append(block.text)
            elif isinstance(block, dict) and 'text' in block:
                parts.append(block['text'])

        return "\n\n".join(parts) if parts else "Empty response"

    def _get_usage_info(self, response: anthropic.types.Message) -> Optional[UsageInfo]:
        """Extract usage information from response."""
        if hasattr(response, 'usage') and response.usage:
            return UsageInfo(
                input_tokens=response.usage.input_tokens,
                output_tokens=response.usage.output_tokens,
                total_tokens=response.usage.input_tokens + response.usage.output_tokens
            )
        return None
