"""Google Gemini provider implementation with MCP tool calling support.

Uses Google's standard Gemini API with manual MCP client management.
"""

import os
import asyncio
import sys
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
    """Google Gemini provider with MCP support via manual tool calling."""

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
        """Send message to Gemini with MCP tool calling.

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
        # Run async tool calling in sync context
        return asyncio.run(
            self._send_message_async(
                messages=messages,
                mcp_servers=mcp_servers,
                model=model,
                max_tokens=max_tokens,
                temperature=temperature,
                uploaded_files=uploaded_files
            )
        )

    async def _send_message_async(
        self,
        messages: List[ChatMessage],
        mcp_servers: List[Dict],
        model: str,
        max_tokens: int,
        temperature: float,
        uploaded_files: Optional[Dict]
    ) -> ChatResponse:
        """Async implementation of message sending with tool calling.

        This implements the agentic loop:
        1. Connect to MCP servers and discover tools
        2. Send message to Gemini with tool declarations
        3. If Gemini calls tools, execute them via MCP
        4. Feed results back to Gemini
        5. Repeat until Gemini responds without tool calls
        """
        try:
            # Import MCP client manager
            # Add parent directory to sys.path to enable absolute imports
            parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if parent_dir not in sys.path:
                sys.path.insert(0, parent_dir)

            from utils.mcp_client import MCPClientManager

            # Connect to MCP servers using context manager
            print(f"DEBUG: Connecting to {len(mcp_servers)} MCP servers", file=sys.stderr)

            async with MCPClientManager(mcp_servers) as mcp_manager:
                # Get Gemini-formatted tools
                gemini_tools = mcp_manager.convert_tools_to_gemini_format()
                print(f"DEBUG: Discovered {len(gemini_tools)} total tools", file=sys.stderr)

                # Convert messages to Gemini format
                gemini_messages = self._convert_messages_to_gemini_format(messages)

                # Build system instruction
                system_instruction = self._build_system_instruction(mcp_servers, uploaded_files)

                # Agentic loop: keep calling until no more tool calls
                # Increased to 30 to handle complex multi-step workflows (e.g., multi-omics integration)
                max_iterations = 30
                iteration = 0
                conversation_history = gemini_messages.copy()
                all_tool_calls = []  # Track all tool calls for trace

                while iteration < max_iterations:
                    iteration += 1
                    print(f"DEBUG: Iteration {iteration}", file=sys.stderr)

                    # Call Gemini with current conversation and tools
                    response = await self._call_gemini(
                        model=model,
                        messages=conversation_history,
                        tools=gemini_tools,
                        system_instruction=system_instruction,
                        max_tokens=max_tokens,
                        temperature=temperature
                    )

                    # Check if response contains tool calls
                    tool_calls, original_parts = self._extract_tool_calls(response)

                    # Debug: log response structure
                    print(f"DEBUG: Response candidates: {len(response.candidates) if hasattr(response, 'candidates') else 0}", file=sys.stderr)
                    if hasattr(response, 'candidates') and response.candidates:
                        cand = response.candidates[0]
                        if hasattr(cand, 'content') and hasattr(cand.content, 'parts'):
                            print(f"DEBUG: Response parts: {len(cand.content.parts)}", file=sys.stderr)
                            for i, part in enumerate(cand.content.parts):
                                print(f"DEBUG: Part {i} type: {type(part).__name__}", file=sys.stderr)
                                if hasattr(part, 'function_call'):
                                    print(f"DEBUG: Part {i} has function_call: {part.function_call}", file=sys.stderr)

                    if not tool_calls:
                        # No more tool calls - we're done
                        print(f"DEBUG: No tool calls, finishing after {iteration} iterations", file=sys.stderr)

                        # Extract final response
                        content = self._format_response(response)
                        usage = self._get_usage_info(response)

                        return ChatResponse(
                            content=content,
                            usage=usage,
                            raw_response=response,
                            tool_calls_metadata=all_tool_calls if all_tool_calls else None
                        )

                    # Execute tool calls
                    print(f"DEBUG: Executing {len(tool_calls)} tool calls", file=sys.stderr)
                    tool_results = await self._execute_tool_calls(mcp_manager, tool_calls)

                    # Store tool calls and results for trace
                    for i, tc in enumerate(tool_calls):
                        all_tool_calls.append({
                            "server_name": tc["name"].split("_", 1)[0] if "_" in tc["name"] else "unknown",
                            "tool_name": tc["name"].split("_", 1)[1] if "_" in tc["name"] else tc["name"],
                            "input": tc["args"],
                            "result": tool_results[i]["response"] if i < len(tool_results) else None
                        })

                    # Add tool calls and results to conversation history
                    # CRITICAL: Use original Part objects directly to preserve thought_signature
                    conversation_history.append({
                        "role": "model",
                        "parts": original_parts  # Use Part objects directly, not wrapped
                    })

                    # Add tool results
                    for result in tool_results:
                        conversation_history.append({
                            "role": "function",
                            "parts": [{
                                "functionResponse": {
                                    "name": result["name"],
                                    "response": result["response"]
                                }
                            }]
                        })

                # Max iterations reached
                print(f"WARNING: Max iterations ({max_iterations}) reached", file=sys.stderr)
                return ChatResponse(
                    content="Error: Maximum tool calling iterations reached. Please try a simpler query.",
                    usage=None,
                    raw_response=None
                )

        except Exception as e:
            import traceback
            traceback.print_exc()
            raise Exception(f"Gemini API error: {str(e)}")

    async def _call_gemini(
        self,
        model: str,
        messages: List[Dict],
        tools: List[Dict],
        system_instruction: str,
        max_tokens: int,
        temperature: float
    ) -> Any:
        """Call Gemini API with tools.

        Args:
            model: Model identifier
            messages: Conversation history
            tools: Tool declarations
            system_instruction: System prompt
            max_tokens: Max output tokens
            temperature: Sampling temperature

        Returns:
            Gemini response
        """
        # Convert tools to Gemini format
        gemini_tools = [
            types.Tool(
                function_declarations=[
                    types.FunctionDeclaration(
                        name=tool["name"],
                        description=tool["description"],
                        parameters=tool["parameters"]
                    )
                    for tool in tools
                ]
            )
        ] if tools else None

        # Build content from messages
        contents = []
        for msg in messages:
            contents.append(types.Content(
                role=msg["role"],
                parts=[types.Part(text=msg["content"])] if "content" in msg else msg.get("parts", [])
            ))

        # Configure function calling to encourage tool usage
        tool_config = None
        if gemini_tools:
            tool_config = types.ToolConfig(
                function_calling_config=types.FunctionCallingConfig(
                    mode="AUTO"  # Let Gemini decide when to call tools
                )
            )

        # Call Gemini
        response = self.client.models.generate_content(
            model=model,
            contents=contents,
            config=types.GenerateContentConfig(
                system_instruction=system_instruction,
                max_output_tokens=max_tokens,
                temperature=temperature,
                tools=gemini_tools,
                tool_config=tool_config
            )
        )

        return response

    def _extract_tool_calls(self, response: Any) -> tuple[List[Dict], List[Any]]:
        """Extract tool calls from Gemini response.

        Args:
            response: Gemini API response

        Returns:
            Tuple of (parsed_tool_calls, original_parts)
            - parsed_tool_calls: List with name and args for execution
            - original_parts: Original Part objects to preserve thought_signature
        """
        parsed_calls = []
        original_parts = []

        if not hasattr(response, 'candidates'):
            return parsed_calls, original_parts

        for candidate in response.candidates:
            if not hasattr(candidate, 'content'):
                continue

            for part in candidate.content.parts:
                if hasattr(part, 'function_call') and part.function_call is not None:
                    fc = part.function_call
                    # Make sure fc has a name attribute and it's not None
                    if hasattr(fc, 'name') and fc.name:
                        parsed_calls.append({
                            "name": fc.name,
                            "args": dict(fc.args) if hasattr(fc, 'args') else {}
                        })
                        # Preserve the ENTIRE Part object (includes thought_signature)
                        original_parts.append(part)

        return parsed_calls, original_parts

    async def _execute_tool_calls(self, mcp_manager, tool_calls: List[Dict]) -> List[Dict]:
        """Execute MCP tool calls.

        Args:
            mcp_manager: MCP client manager instance
            tool_calls: List of tool calls from Gemini

        Returns:
            List of tool results
        """
        results = []

        for tool_call in tool_calls:
            # Parse server name and tool name (format: servername_toolname)
            full_name = tool_call["name"]
            parts = full_name.split("_", 1)

            if len(parts) != 2:
                print(f"ERROR: Invalid tool name format: {full_name}", file=sys.stderr)
                results.append({
                    "name": full_name,
                    "response": {"error": f"Invalid tool name format: {full_name}"}
                })
                continue

            server_name, tool_name = parts

            # Call the tool via MCP manager
            try:
                result = await mcp_manager.call_tool(
                    server_name=server_name,
                    tool_name=tool_name,
                    arguments=tool_call["args"]
                )

                # Format result for Gemini
                results.append({
                    "name": full_name,
                    "response": {
                        "content": result["content"],
                        "isError": result.get("isError", False)
                    }
                })

            except Exception as e:
                print(f"ERROR: Tool execution failed: {e}", file=sys.stderr)
                results.append({
                    "name": full_name,
                    "response": {"error": str(e)}
                })

        return results

    def get_provider_name(self) -> str:
        """Get provider name."""
        return "Gemini"

    def get_model_display_name(self, model: str) -> str:
        """Get human-readable model name."""
        model_names = {
            "gemini-3-flash-preview": "Gemini 3 Flash (Preview)",
            "gemini-3-pro-preview": "Gemini 3 Pro (Preview)",
            "gemini-2.5-flash": "Gemini 2.5 Flash",
            "gemini-2.5-flash-preview-09-2025": "Gemini 2.5 Flash (Preview)",
        }
        return model_names.get(model, model)

    def format_mcp_servers(self, server_configs: List[Dict]) -> List[Dict]:
        """Format MCP servers for Gemini (pass-through).

        Note: With the SSE-based approach, MCP server formatting
        is handled by the MCPClientManager, so this just returns
        the configs as-is.

        Args:
            server_configs: Standard MCP server configs

        Returns:
            Server configs (unchanged)
        """
        return server_configs

    def is_available(self) -> bool:
        """Check if Gemini provider is available."""
        return GEMINI_AVAILABLE and bool(self.api_key)

    def _convert_messages_to_gemini_format(self, messages: List[ChatMessage]) -> List[Dict]:
        """Convert ChatMessage objects to Gemini format.

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
                "content": msg.content
            })

        return gemini_messages

    def _build_system_instruction(
        self,
        mcp_servers: List[Dict],
        uploaded_files: Optional[Dict] = None
    ) -> str:
        """Build system instruction for Gemini."""
        server_descriptions = [
            f"- {server['name']}: Available for bioinformatics analysis"
            for server in mcp_servers
        ]

        system_instruction = f"""You are a bioinformatics assistant with access to specialized analysis tools via MCP servers.

Available MCP servers:
{chr(10).join(server_descriptions)}

When to use tools:
- For ANY analysis request (pathway enrichment, spatial analysis, genomic analysis, etc.) - call the appropriate tool
- For validating or inspecting data files - use the validation tools
- For computational tasks - execute the relevant tool
- Prefer calling tools over providing theoretical responses

When NOT to use tools:
- Simple questions about capabilities or available tools
- Clarifying questions about parameters or requirements
- Explaining tool outputs you've already received

After receiving tool results:
- Interpret and explain the results clearly to the user
- Provide actionable insights based on the analysis output
- Stop calling tools once you have sufficient information to answer the query"""

        # Add file information if present
        if uploaded_files and len(uploaded_files) > 0:
            file_count = len(uploaded_files)
            file_list = []
            for filename, file_info in uploaded_files.items():
                file_name = file_info.get('original_name', filename)
                file_path = file_info.get('path', filename)
                source = file_info.get('source', 'local')

                # For GCS files, show the full URI
                if source == 'gcs':
                    file_list.append(f"- {file_name} (GCS URI: {file_path})")
                else:
                    file_list.append(f"- {file_name} (local path: {file_path})")

            system_instruction += f"""

UPLOADED FILES:
The user has uploaded {file_count} file(s):
{chr(10).join(file_list)}

IMPORTANT: When calling MCP tools with these files:
- For GCS files: Pass the full GCS URI (e.g., gs://bucket/path/file.csv) as the file parameter
- For local files: Pass the local file path as the file parameter
- The MCP servers running on Cloud Run can access GCS URIs directly"""

        return system_instruction

    def _format_response(self, response) -> str:
        """Format Gemini response for display.

        Args:
            response: Gemini response object

        Returns:
            Formatted text response
        """
        if not hasattr(response, 'candidates'):
            return str(response)

        parts = []
        for candidate in response.candidates:
            if hasattr(candidate, 'content'):
                for part in candidate.content.parts:
                    if hasattr(part, 'text'):
                        parts.append(part.text)

        return "\n\n".join(parts) if parts else "No response generated"

    def _get_usage_info(self, response) -> Optional[UsageInfo]:
        """Extract usage information from Gemini response.

        Args:
            response: Gemini response object

        Returns:
            UsageInfo if available, None otherwise
        """
        if hasattr(response, 'usage_metadata'):
            usage = response.usage_metadata
            return UsageInfo(
                input_tokens=getattr(usage, 'prompt_token_count', 0),
                output_tokens=getattr(usage, 'candidates_token_count', 0),
                total_tokens=getattr(usage, 'total_token_count', 0)
            )

        return None
