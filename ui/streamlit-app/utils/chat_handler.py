"""Chat handler for Claude API with MCP servers.

Handles message formatting, API calls, and response parsing.
"""

import os
from typing import List, Dict, Optional
import anthropic


class ChatHandler:
    """Handle chat interactions with Claude API and MCP servers."""

    def __init__(self, api_key: Optional[str] = None):
        """Initialize chat handler.

        Args:
            api_key: Anthropic API key (defaults to ANTHROPIC_API_KEY env var)
        """
        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError("ANTHROPIC_API_KEY not found in environment")

        self.client = anthropic.Anthropic(api_key=self.api_key)

    def send_message(
        self,
        messages: List[Dict[str, str]],
        mcp_servers: List[Dict],
        model: str = "claude-sonnet-4-5",
        max_tokens: int = 4096,
        temperature: float = 1.0,
        uploaded_files: Optional[Dict] = None
    ) -> anthropic.types.Message:
        """Send message to Claude with MCP servers enabled.

        Args:
            messages: Chat messages history
            mcp_servers: MCP server configurations
            model: Claude model to use
            max_tokens: Maximum response tokens
            temperature: Sampling temperature
            uploaded_files: Dict of uploaded files with paths and metadata

        Returns:
            Claude API response message
        """
        # Get tools config for selected servers
        tools = [
            {"type": "mcp_toolset", "mcp_server_name": server["name"]}
            for server in mcp_servers
        ]

        # Build system prompt instructing Claude to use MCP tools
        server_descriptions = []
        for server in mcp_servers:
            server_descriptions.append(f"- {server['name']}: Available for bioinformatics analysis")

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

        # Add uploaded file information to system prompt
        if uploaded_files and len(uploaded_files) > 0:
            file_descriptions = []
            for filename, file_info in uploaded_files.items():
                metadata = file_info['metadata']
                path = file_info['path']
                original_name = file_info.get('original_name', filename)

                # Check if this is a GCS file or local file
                is_gcs = metadata.get('source') == 'gcs'

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

                    # For small text files, include content inline
                    content = None

                    # Check if GCS content was pre-loaded
                    if '_gcs_content' in metadata:
                        content = metadata['_gcs_content']
                    # For local files under 50KB, read content
                    elif not is_gcs and metadata.get('size_bytes', 0) < 50000:
                        try:
                            with open(path, 'r') as f:
                                content = f.read()
                        except Exception:
                            pass  # Skip if can't read

                    # Include content if available
                    if content:
                        # Truncate if too long
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

        try:
            # Use beta.messages API for MCP support
            # MCP is built into beta API - no additional beta header needed
            response = self.client.beta.messages.create(
                model=model,
                max_tokens=max_tokens,
                temperature=temperature,
                messages=messages,
                mcp_servers=mcp_servers,
                tools=tools,
                system=system_prompt
            )
            return response

        except anthropic.APIError as e:
            raise Exception(f"Claude API error: {str(e)}")
        except Exception as e:
            raise Exception(f"Unexpected error: {str(e)}")

    def format_response(self, response: anthropic.types.Message) -> str:
        """Format Claude response for display.

        Args:
            response: Claude API response

        Returns:
            Formatted response text
        """
        if not response.content:
            return "No response from Claude"

        # Extract text from content blocks
        parts = []
        for block in response.content:
            if hasattr(block, 'text'):
                parts.append(block.text)
            elif isinstance(block, dict) and 'text' in block:
                parts.append(block['text'])

        return "\n\n".join(parts) if parts else "Empty response"

    def get_usage_info(self, response: anthropic.types.Message) -> Dict:
        """Extract usage information from response.

        Args:
            response: Claude API response

        Returns:
            Dict with token usage info
        """
        usage = response.usage if hasattr(response, 'usage') else None
        if usage:
            return {
                "input_tokens": usage.input_tokens,
                "output_tokens": usage.output_tokens,
                "total_tokens": usage.input_tokens + usage.output_tokens
            }
        return {}
