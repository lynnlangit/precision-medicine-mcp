# LLM Provider Architecture

This directory contains the provider abstraction layer that enables multiple LLM providers (Claude, Gemini) to work with the same MCP server infrastructure.

## Overview

The provider system uses an abstract base class pattern to define a common interface for different LLM providers, each with their own MCP integration approach.

```
┌─────────────────────────────────────────┐
│         LLMProvider (Abstract)          │
│  - send_message()                       │
│  - get_provider_name()                  │
│  - format_mcp_servers()                 │
└─────────────────────────────────────────┘
           ▲                    ▲
           │                    │
    ┌──────┴──────┐      ┌─────┴──────┐
    │   Claude    │      │   Gemini   │
    │  Provider   │      │  Provider  │
    └─────────────┘      └────────────┘
```

## Files

### `base.py` - Abstract Base Class

Defines the interface all providers must implement:

**Data Classes:**
- `ChatMessage` - Standardized message format
- `ChatResponse` - Standardized response format
- `UsageInfo` - Token usage tracking

**Abstract Methods:**
- `send_message()` - Send chat message with MCP servers
- `get_provider_name()` - Return provider name for display
- `get_model_display_name()` - Format model names
- `format_mcp_servers()` - Convert MCP configs to provider format

### `anthropic_provider.py` - Claude Integration

Uses Anthropic's native MCP support via the Claude API.

**Architecture:**
```
AnthropicProvider → Claude API (with mcp_servers) → Response
```

**Key Features:**
- Native MCP server integration
- Claude API handles all orchestration
- Direct tool discovery and calling
- No custom client needed

**Implementation:**
```python
response = self.client.messages.create(
    model=model,
    messages=messages,
    mcp_servers=formatted_servers,  # Native support
    max_tokens=max_tokens
)
```

### `gemini_provider.py` - Gemini Integration

Custom SSE-based MCP client with manual orchestration.

**Architecture:**
```
GeminiProvider → MCP Client Manager → SSE Connections → Cloud Run Servers
              ↓                                               ↓
         Gemini API ← Tool Results ← ← ← ← Tool Execution ← ←
```

**Key Features:**
- Direct SSE connections to MCP servers
- Agentic tool calling loop
- Schema cleaning for Gemini compatibility
- Google Cloud authentication
- Preserves thought signatures for multi-turn calls

**Implementation Flow:**
1. **Connect** - Establish SSE connections to all MCP servers
2. **Discover** - Fetch available tools from each server
3. **Clean** - Convert MCP schemas to Gemini format
4. **Loop** - Agentic tool calling:
   - Gemini decides which tools to call
   - Execute tools via MCP client
   - Feed results back to Gemini
   - Repeat until complete
5. **Cleanup** - Close all connections

**Why Custom Client:**
- Gemini's Interactions API doesn't support remote MCP servers
- No built-in tool configuration or forcing
- Requires manual tool orchestration
- Needs direct control over connection lifecycle

### `__init__.py` - Provider Factory

Provides helper functions for provider discovery and instantiation.

**Functions:**
- `get_provider(name, api_key)` - Create provider instance
- `get_available_providers()` - List providers with status

**Usage:**
```python
from providers import get_provider

# Create provider
provider = get_provider("gemini", api_key="your_key")

# Send message
response = provider.send_message(
    messages=chat_messages,
    mcp_servers=mcp_configs,
    model="gemini-3-flash-preview"
)
```

## Adding a New Provider

To add support for a new LLM provider:

1. **Create provider file**: `providers/newprovider_provider.py`

2. **Inherit from base class**:
```python
from .base import LLMProvider, ChatMessage, ChatResponse

class NewProvider(LLMProvider):
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.client = NewProviderClient(api_key)
```

3. **Implement required methods**:
```python
def send_message(self, messages, mcp_servers, model, ...):
    # Your implementation
    pass

def get_provider_name(self) -> str:
    return "NewProvider"

def format_mcp_servers(self, server_configs):
    # Convert to provider-specific format
    pass
```

4. **Register in factory** (`__init__.py`):
```python
from .newprovider_provider import NewProvider

def get_provider(provider_name: str, api_key: Optional[str] = None):
    if provider_name == "newprovider":
        return NewProvider(api_key=api_key)
```

5. **Update available providers**:
```python
def get_available_providers() -> dict:
    providers["newprovider"] = {
        "name": "NewProvider",
        "available": bool(os.getenv("NEWPROVIDER_API_KEY")),
        "models": [{"id": "model-1", "name": "Model 1"}]
    }
```

## MCP Integration Patterns

### Pattern 1: Native MCP Support (Claude)

**When to use:**
- Provider has built-in MCP support
- API handles tool orchestration automatically

**Pros:**
- Simple implementation
- Provider handles all complexity
- No custom client needed

**Example:**
```python
response = self.client.messages.create(
    mcp_servers=formatted_servers  # Native parameter
)
```

### Pattern 2: SSE Client + Manual Orchestration (Gemini)

**When to use:**
- Provider lacks MCP support
- Need direct control over tool calling
- Provider requires custom schema format

**Pros:**
- Full control over tool execution
- Can adapt any MCP server
- Custom error handling

**Cons:**
- More complex implementation
- Manual connection management
- Custom schema conversion needed

**Example:**
```python
async with MCPClientManager(mcp_servers) as mcp_manager:
    tools = mcp_manager.convert_tools_to_gemini_format()
    # Manual orchestration loop
    while has_tool_calls:
        results = await mcp_manager.call_tool(...)
```

### Pattern 3: HTTP API Bridge (Future)

**When to use:**
- Provider requires REST/HTTP tool format
- No SSE support

**Implementation:**
- Convert MCP SSE to HTTP requests
- Polling or webhooks for results
- Request/response transformation

## Schema Compatibility

Different providers have different requirements for function calling schemas:

**Gemini Restrictions:**
- No `additionalProperties`, `anyOf`, `oneOf`, `allOf`
- No `$ref`, `$schema`, `$id`
- Requires `thought_signature` in multi-turn calls
- Must preserve complete Part objects

**Claude:**
- Full JSON Schema support
- Native MCP tool format
- No special requirements

**Solution:** Provider-specific schema cleaning in `format_mcp_servers()` or dedicated cleaning functions.

## Testing

Test each provider with:

```python
# Test basic functionality
provider = get_provider("gemini")
response = provider.send_message(
    messages=[ChatMessage(role="user", content="test")],
    mcp_servers=[{"name": "test", "url": "..."}],
    model="gemini-3-flash-preview"
)

# Verify response format
assert isinstance(response, ChatResponse)
assert response.content
```

## Future Enhancements

- [ ] OpenAI provider (HTTP-based MCP adapter)
- [ ] Azure OpenAI provider
- [ ] AWS Bedrock provider
- [ ] Provider-specific error handling
- [ ] Retry logic with exponential backoff
- [ ] Rate limiting per provider
- [ ] Provider health checks
- [ ] Cost estimation per provider
