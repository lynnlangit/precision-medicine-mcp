"""Streamlit Chat Interface for MCP Servers

A visual interface for testing deployed MCP servers on GCP Cloud Run.
Provides a Claude Desktop-like experience for bioinformatics workflows.
"""

import streamlit as st
import os
import uuid
import tempfile
from typing import List, Dict
from datetime import datetime
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

from utils import (
    MCP_SERVERS,
    get_server_config,
    get_tools_config,
    get_server_categories,
    EXAMPLE_PROMPTS,
    ChatHandler,
    build_orchestration_trace,
    render_trace,
    render_trace_export
)

# Import provider system
from providers import (
    get_provider,
    get_available_providers,
    ChatMessage,
    GEMINI_AVAILABLE
)

# Import authentication and audit logging
from utils.auth import require_authentication, display_user_info, display_logout_button
from utils.audit_logger import get_audit_logger

# Import file validation and GCS handling
from utils.file_validator import validate_uploaded_file, sanitize_filename
from utils.gcs_handler import (
    is_gcs_path,
    validate_gcs_uri,
    get_gcs_file_metadata,
    get_gcs_file_content
)

# Page configuration
st.set_page_config(
    page_title="MCP Chat - Precision Medicine",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .stChatMessage {
        padding: 1rem;
        border-radius: 0.5rem;
    }
    .server-card {
        padding: 1rem;
        border-radius: 0.5rem;
        border: 1px solid #ddd;
        margin: 0.5rem 0;
    }
    .production {
        border-left: 4px solid #28a745;
    }
    .mock {
        border-left: 4px solid #ffc107;
    }
</style>
""", unsafe_allow_html=True)


def is_cloud_run() -> bool:
    """Detect if running on GCP Cloud Run.

    Returns:
        bool: True if deployed on Cloud Run, False for local development
    """
    return os.getenv("K_SERVICE") is not None


def initialize_session_state():
    """Initialize Streamlit session state."""
    if "messages" not in st.session_state:
        st.session_state.messages = []
    if "selected_servers" not in st.session_state:
        st.session_state.selected_servers = ["spatialtools", "multiomics", "fgbio"]

    # Initialize provider selection (Cloud Run only)
    if "llm_provider" not in st.session_state:
        st.session_state.llm_provider = "claude"  # Default to Claude

    if "llm_model" not in st.session_state:
        st.session_state.llm_model = "claude-sonnet-4-20250514"  # Default model

    # Initialize provider instance
    if "provider_instance" not in st.session_state:
        try:
            st.session_state.provider_instance = get_provider(
                provider_name=st.session_state.llm_provider
            )
        except (ValueError, ImportError) as e:
            st.session_state.provider_instance = None
            if is_cloud_run():
                st.error(f"Failed to initialize {st.session_state.llm_provider} provider: {str(e)}")

    # Keep chat_handler for backward compatibility (local development)
    if "chat_handler" not in st.session_state:
        api_key = os.getenv("ANTHROPIC_API_KEY")
        if api_key:
            try:
                st.session_state.chat_handler = ChatHandler(api_key)
            except ValueError:
                st.session_state.chat_handler = None
        else:
            st.session_state.chat_handler = None

    # Initialize session tracking
    if "session_id" not in st.session_state:
        st.session_state.session_id = f"sess_{uuid.uuid4().hex[:12]}"
    if "session_start_time" not in st.session_state:
        st.session_state.session_start_time = datetime.utcnow()
    if "total_queries" not in st.session_state:
        st.session_state.total_queries = 0
    if "total_tokens" not in st.session_state:
        st.session_state.total_tokens = 0
    if "total_cost" not in st.session_state:
        st.session_state.total_cost = 0.0

    # Initialize audit logger
    if "audit_logger" not in st.session_state:
        st.session_state.audit_logger = get_audit_logger()

    # Initialize trace storage
    if "traces" not in st.session_state:
        st.session_state.traces = {}  # message_index -> OrchestrationTrace

    # Initialize uploaded files storage
    if "uploaded_files" not in st.session_state:
        st.session_state.uploaded_files = {}  # filename -> {'path': temp_path, 'metadata': dict}


def prepare_file_for_mcp(uploaded_file, metadata: Dict) -> str:
    """Prepare uploaded file for MCP server access.

    Writes validated file to a temporary directory and returns the path
    for MCP servers to access.

    Args:
        uploaded_file: Streamlit UploadedFile object
        metadata: File metadata dict from validation

    Returns:
        str: Absolute path to the temporary file

    Example:
        >>> temp_path = prepare_file_for_mcp(uploaded_file, metadata)
        >>> # Pass temp_path to MCP tool call
    """
    # Create a temp directory if it doesn't exist
    temp_dir = Path(tempfile.gettempdir()) / "mcp_uploads"
    temp_dir.mkdir(exist_ok=True)

    # Use sanitized filename
    sanitized_name = metadata['sanitized_filename']
    temp_path = temp_dir / sanitized_name

    # Write file content
    content = uploaded_file.read()
    uploaded_file.seek(0)  # Reset file pointer

    with open(temp_path, 'wb') as f:
        f.write(content)

    return str(temp_path)


def render_sidebar():
    """Render sidebar with server selection and settings."""
    with st.sidebar:
        st.title("🧬 MCP Chat")
        st.markdown("---")

        # User info (SSO authentication)
        if "user" in st.session_state:
            display_user_info(st.session_state.user, location="sidebar")
            display_logout_button()

        st.markdown("---")

        # API Key status
        if st.session_state.chat_handler:
            st.success("✅ API Key Configured")
        else:
            st.error("❌ API Key Missing")
            st.info("Set ANTHROPIC_API_KEY environment variable")

        st.markdown("---")

        # Provider selection (Cloud Run only)
        if is_cloud_run():
            st.subheader("LLM Provider")

            available_providers = get_available_providers()

            # Build provider options
            provider_options = []
            provider_mapping = {}
            for provider_key, provider_info in available_providers.items():
                if provider_info["available"]:
                    display_name = provider_info["name"]
                    provider_options.append(display_name)
                    provider_mapping[display_name] = provider_key

            if not provider_options:
                st.error("No LLM providers available. Configure API keys.")
                st.stop()

            # Provider dropdown
            current_provider_display = available_providers[st.session_state.llm_provider]["name"]
            selected_provider_display = st.selectbox(
                "Select Provider",
                options=provider_options,
                index=provider_options.index(current_provider_display) if current_provider_display in provider_options else 0,
                help="Choose LLM provider for analysis"
            )

            # Update provider if changed
            selected_provider_key = provider_mapping[selected_provider_display]
            if selected_provider_key != st.session_state.llm_provider:
                st.session_state.llm_provider = selected_provider_key
                try:
                    st.session_state.provider_instance = get_provider(provider_name=selected_provider_key)
                    st.success(f"Switched to {selected_provider_display}")
                    st.rerun()
                except (ValueError, ImportError) as e:
                    st.error(f"Failed to initialize {selected_provider_display}: {str(e)}")
                    st.stop()

            st.markdown("---")

        # Model selection
        st.subheader("Model Settings")

        # Dynamic model selection based on provider (Cloud Run)
        if is_cloud_run() and st.session_state.provider_instance:
            available_providers = get_available_providers()
            provider_info = available_providers.get(st.session_state.llm_provider, {})
            available_models = provider_info.get("models", [])

            if available_models:
                model_options = [m["id"] for m in available_models]
                model_labels = [m["name"] for m in available_models]

                # Create mapping for display
                model_display_map = {m["name"]: m["id"] for m in available_models}

                current_model = st.session_state.llm_model
                try:
                    current_index = model_options.index(current_model)
                except ValueError:
                    current_index = 0
                    st.session_state.llm_model = model_options[0]

                selected_model_label = st.selectbox(
                    f"{provider_info['name']} Model",
                    options=model_labels,
                    index=current_index,
                    help=f"Select {provider_info['name']} model"
                )
                model = model_display_map[selected_model_label]
                st.session_state.llm_model = model
            else:
                model = provider_info.get("default_model", "claude-sonnet-4-20250514")
                st.session_state.llm_model = model
        else:
            # Local development - Claude only
            model = st.selectbox(
                "Claude Model",
                options=["claude-sonnet-4-5", "claude-opus-4-5", "claude-haiku-4"],
                index=0,
                help="Select Claude model. Sonnet recommended for most tasks."
            )

        max_tokens = st.slider(
            "Max Tokens",
            min_value=1024,
            max_value=8192,
            value=4096,
            step=512,
            help="Maximum response length"
        )

        st.markdown("---")

        # MCP Server Selection
        st.subheader("MCP Servers")

        categories = get_server_categories()

        # Production servers
        st.markdown("**Production Servers** (Real Analysis)")
        for server_name in categories["Production Servers (Real Analysis)"]:
            server = MCP_SERVERS[server_name]
            selected = st.checkbox(
                f"{server_name} ({server['tools_count']} tools)",
                value=server_name in st.session_state.selected_servers,
                key=f"checkbox_{server_name}",
                help=server['description']
            )
            if selected and server_name not in st.session_state.selected_servers:
                st.session_state.selected_servers.append(server_name)
            elif not selected and server_name in st.session_state.selected_servers:
                st.session_state.selected_servers.remove(server_name)

        st.markdown("**Mock Servers** (Demo Only)")
        with st.expander("Show Mock Servers"):
            for server_name in categories["Mock Servers (Workflow Demo)"]:
                server = MCP_SERVERS[server_name]
                selected = st.checkbox(
                    f"{server_name} ({server['tools_count']} tools)",
                    value=server_name in st.session_state.selected_servers,
                    key=f"checkbox_{server_name}",
                    help=server['description']
                )
                if selected and server_name not in st.session_state.selected_servers:
                    st.session_state.selected_servers.append(server_name)
                elif not selected and server_name in st.session_state.selected_servers:
                    st.session_state.selected_servers.remove(server_name)

        # Show selected count
        st.info(f"**{len(st.session_state.selected_servers)}** servers selected")

        st.markdown("---")

        # Orchestration Trace Settings
        st.subheader("🔍 Orchestration Trace")
        show_trace = st.toggle(
            "Show trace for responses",
            value=False,
            help="Display which MCP servers were called and in what order"
        )

        trace_style = "log"  # Default
        if show_trace:
            trace_style = st.selectbox(
                "Trace style",
                options=["log", "cards", "timeline", "mermaid"],
                format_func=lambda x: {
                    "log": "📝 Log View",
                    "cards": "🎴 Card View",
                    "timeline": "📈 Timeline View",
                    "mermaid": "📊 Sequence Diagram"
                }.get(x, x),
                help="Choose how to display the orchestration trace"
            )

        st.markdown("---")

        # File Upload Section
        st.subheader("📁 File Upload")
        st.caption("Upload bioinformatics data files")

        uploaded_files = st.file_uploader(
            "Choose files",
            type=['fasta', 'fa', 'fna', 'fastq', 'fq', 'vcf', 'gff', 'gtf', 'bed',
                  'csv', 'tsv', 'tab', 'txt', 'json', 'h5ad', 'h5',
                  'png', 'jpg', 'jpeg', 'tiff', 'tif'],
            accept_multiple_files=True,
            help="Upload sequence files, annotations, tabular data, or images"
        )

        if uploaded_files:
            st.caption(f"{len(uploaded_files)} file(s) selected")

            for uploaded_file in uploaded_files:
                # Validate each file
                is_valid, errors, metadata = validate_uploaded_file(uploaded_file)

                if is_valid:
                    # Store validated file
                    temp_path = prepare_file_for_mcp(uploaded_file, metadata)
                    st.session_state.uploaded_files[metadata['sanitized_filename']] = {
                        'path': temp_path,
                        'metadata': metadata,
                        'original_name': uploaded_file.name
                    }

                    # Show success with file info
                    with st.expander(f"✅ {uploaded_file.name}", expanded=False):
                        st.success("Valid file")
                        col1, col2 = st.columns(2)
                        col1.metric("Size", f"{metadata['size_mb']:.2f} MB")
                        col2.metric("Type", metadata['extension'])
                        if not metadata.get('is_binary', False):
                            st.caption(f"Lines: {metadata.get('line_count', 'N/A')}")
                else:
                    # Show validation errors
                    with st.expander(f"❌ {uploaded_file.name}", expanded=True):
                        st.error("Invalid file")
                        for error in errors:
                            st.write(f"- {error}")

        # GCS Path Input Section
        st.caption("Or provide GCS bucket path")

        gcs_uri = st.text_input(
            "GCS URI",
            placeholder="gs://bucket-name/path/to/file.fastq",
            help="Enter a Google Cloud Storage URI to access files from GCS",
            key="gcs_uri_input"
        )

        if gcs_uri and gcs_uri.strip():
            # Validate GCS URI
            is_valid, error = validate_gcs_uri(gcs_uri)

            if is_valid:
                # Get file metadata
                success, metadata, meta_error = get_gcs_file_metadata(gcs_uri, use_mock=True)

                if success:
                    # Store GCS file reference
                    st.session_state.uploaded_files[metadata['sanitized_filename']] = {
                        'path': gcs_uri,  # Store GCS path directly
                        'metadata': metadata,
                        'original_name': metadata['filename'],
                        'source': 'gcs'
                    }

                    # Show success
                    with st.expander(f"✅ {metadata['filename']} (GCS)", expanded=False):
                        st.success("Valid GCS URI")
                        col1, col2 = st.columns(2)
                        col1.metric("Bucket", metadata['bucket'])
                        col2.metric("Type", metadata['extension'])
                        st.caption(f"Path: {metadata['blob_path']}")

                        # Try to get content for small files
                        if not metadata.get('is_binary', False):
                            success_content, content, content_error = get_gcs_file_content(gcs_uri, max_size_bytes=50000)
                            if success_content:
                                st.caption(f"✅ Content loaded ({len(content)} chars)")
                                # Store content in metadata for inline inclusion
                                metadata['_gcs_content'] = content
                            else:
                                st.caption(f"⚠️ Content not loaded: {content_error}")
                else:
                    st.error(f"Error: {meta_error}")
            else:
                st.error(f"Invalid GCS URI: {error}")

        # Show currently uploaded files
        if st.session_state.uploaded_files:
            st.caption(f"📎 {len(st.session_state.uploaded_files)} file(s) available for MCP")

            # Button to clear all uploaded files
            if st.button("Clear Files", key="clear_files"):
                # Clean up temp files
                for file_info in st.session_state.uploaded_files.values():
                    try:
                        if os.path.exists(file_info['path']):
                            os.remove(file_info['path'])
                    except Exception:
                        pass
                st.session_state.uploaded_files = {}
                st.rerun()

        st.markdown("---")

        # Example prompts
        st.subheader("Example Prompts")
        selected_example = st.selectbox(
            "Load an example:",
            options=[""] + list(EXAMPLE_PROMPTS.keys()),
            index=0
        )

        if selected_example and st.button("Load Prompt"):
            st.session_state.example_prompt = EXAMPLE_PROMPTS[selected_example]

        st.markdown("---")

        # Session statistics
        st.subheader("Session Stats")
        col1, col2 = st.columns(2)
        col1.metric("Queries", st.session_state.total_queries)
        col2.metric("Tokens", f"{st.session_state.total_tokens:,}")
        st.metric("Cost", f"${st.session_state.total_cost:.4f}")

        # Session duration
        if st.session_state.session_start_time:
            duration = datetime.utcnow() - st.session_state.session_start_time
            duration_mins = int(duration.total_seconds() / 60)
            st.caption(f"Session: {duration_mins} min")

        st.markdown("---")

        # Clear chat button
        if st.button("🗑️ Clear Chat", use_container_width=True):
            st.session_state.messages = []
            st.session_state.traces = {}  # Clear traces too
            st.rerun()

    return model, max_tokens, show_trace, trace_style


def render_server_status():
    """Render server status cards."""
    st.subheader("Active MCP Servers")

    if not st.session_state.selected_servers:
        st.warning("⚠️ No servers selected. Select servers from the sidebar.")
        return

    cols = st.columns(3)
    for idx, server_name in enumerate(st.session_state.selected_servers):
        server = MCP_SERVERS[server_name]
        with cols[idx % 3]:
            status_class = "production" if server["status"] == "production" else "mock"
            st.markdown(f"""
            <div class="server-card {status_class}">
                <h4>{server_name}</h4>
                <p>{server['description']}</p>
                <small>{server['tools_count']} tools • {server['status']}</small>
            </div>
            """, unsafe_allow_html=True)


def render_chat_history(show_trace: bool = False, trace_style: str = "log"):
    """Render chat message history.

    Args:
        show_trace: Whether to show orchestration traces
        trace_style: Style for rendering traces
    """
    for i, message in enumerate(st.session_state.messages):
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

            # Show usage info if available
            if "usage" in message and message["usage"]:
                with st.expander("Token Usage"):
                    usage = message["usage"]
                    col1, col2, col3 = st.columns(3)
                    col1.metric("Input", usage.get("input_tokens", 0))
                    col2.metric("Output", usage.get("output_tokens", 0))
                    col3.metric("Total", usage.get("total_tokens", 0))

            # Show trace for assistant messages if enabled
            if message["role"] == "assistant" and show_trace:
                if i in st.session_state.traces:
                    trace = st.session_state.traces[i]
                    render_trace(trace, style=trace_style)

                    # Add export buttons if trace has data
                    if trace.tool_calls:
                        render_trace_export(trace, trace_id=i)


def handle_user_input(prompt: str, model: str, max_tokens: int):
    """Handle user input and get response from Claude.

    Args:
        prompt: User's message
        model: Claude model to use
        max_tokens: Maximum tokens for response
    """
    # Get user info for audit logging
    user = st.session_state.get("user")
    if not user:
        # Development mode - create mock user
        user = {
            "email": "developer@localhost",
            "user_id": "dev_user",
            "display_name": "Developer"
        }

    # Add user message to history
    st.session_state.messages.append({"role": "user", "content": prompt})

    # Display user message
    with st.chat_message("user"):
        st.markdown(prompt)

    # Get MCP server config
    mcp_servers = get_server_config(st.session_state.selected_servers)

    if not mcp_servers:
        with st.chat_message("assistant"):
            error_msg = "⚠️ No MCP servers selected. Please select at least one server from the sidebar."
            st.error(error_msg)
            st.session_state.messages.append({"role": "assistant", "content": error_msg})
        return

    # Log the query (audit trail)
    st.session_state.audit_logger.log_mcp_query(
        user_email=user["email"],
        user_id=user["user_id"],
        servers=st.session_state.selected_servers,
        prompt=prompt,
        model=model,
        session_id=st.session_state.session_id
    )

    # Track query start time
    import time
    query_start_time = datetime.utcnow()
    start_time = time.time()

    # Show thinking indicator
    with st.chat_message("assistant"):
        with st.spinner("Thinking..."):
            try:
                # Use provider abstraction (Cloud Run) or ChatHandler (local)
                if is_cloud_run() and st.session_state.provider_instance:
                    # Convert messages to ChatMessage objects
                    chat_messages = [
                        ChatMessage(role=msg["role"], content=msg["content"])
                        for msg in st.session_state.messages
                    ]

                    # Send via provider abstraction
                    chat_response = st.session_state.provider_instance.send_message(
                        messages=chat_messages,
                        mcp_servers=mcp_servers,
                        model=model,
                        max_tokens=max_tokens,
                        uploaded_files=st.session_state.uploaded_files
                    )

                    # Calculate response time
                    duration_ms = (time.time() - start_time) * 1000

                    # Extract response content and usage
                    response_text = chat_response.content
                    usage = None
                    if chat_response.usage:
                        usage = {
                            "input_tokens": chat_response.usage.input_tokens,
                            "output_tokens": chat_response.usage.output_tokens,
                            "total_tokens": chat_response.usage.total_tokens
                        }

                    # Keep raw response for trace building
                    response = chat_response.raw_response

                else:
                    # Local development - use existing ChatHandler
                    response = st.session_state.chat_handler.send_message(
                        messages=[{"role": msg["role"], "content": msg["content"]}
                                 for msg in st.session_state.messages],
                        mcp_servers=mcp_servers,
                        model=model,
                        max_tokens=max_tokens,
                        uploaded_files=st.session_state.uploaded_files
                    )

                    # Calculate response time
                    duration_ms = (time.time() - start_time) * 1000

                    # Format response (existing logic)
                    response_text = st.session_state.chat_handler.format_response(response)
                    usage = st.session_state.chat_handler.get_usage_info(response)

                # Calculate query duration
                query_duration = (datetime.utcnow() - query_start_time).total_seconds()

                # Update session statistics
                st.session_state.total_queries += 1
                if usage:
                    st.session_state.total_tokens += usage.get("total_tokens", 0)
                    # Estimate cost (Sonnet pricing)
                    input_cost = usage.get("input_tokens", 0) * 0.003 / 1000
                    output_cost = usage.get("output_tokens", 0) * 0.015 / 1000
                    estimated_cost = input_cost + output_cost
                    st.session_state.total_cost += estimated_cost

                    # Log the response (audit trail)
                    st.session_state.audit_logger.log_mcp_response(
                        user_email=user["email"],
                        user_id=user["user_id"],
                        servers=st.session_state.selected_servers,
                        response_length=len(response_text),
                        input_tokens=usage.get("input_tokens", 0),
                        output_tokens=usage.get("output_tokens", 0),
                        total_tokens=usage.get("total_tokens", 0),
                        estimated_cost=estimated_cost,
                        duration_seconds=query_duration,
                        session_id=st.session_state.session_id
                    )

                # Display response
                st.markdown(response_text)

                # Show usage info
                if usage:
                    with st.expander("Token Usage"):
                        col1, col2, col3 = st.columns(3)
                        col1.metric("Input", usage.get("input_tokens", 0))
                        col2.metric("Output", usage.get("output_tokens", 0))
                        col3.metric("Total", usage.get("total_tokens", 0))

                        # Show estimated cost
                        if usage.get("total_tokens", 0) > 0:
                            st.caption(f"Est. cost: ${estimated_cost:.4f}")

                # Build orchestration trace
                message_index = len(st.session_state.messages)  # Index where this message will be stored
                trace = build_orchestration_trace(
                    query=prompt,
                    response=response,
                    messages=[{"role": msg["role"], "content": msg["content"]}
                             for msg in st.session_state.messages],
                    duration_ms=duration_ms,
                    tokens=usage.get("total_tokens", 0),
                    cost_usd=estimated_cost if usage else 0
                )

                # Store trace in session state
                st.session_state.traces[message_index] = trace

                # Add to history
                st.session_state.messages.append({
                    "role": "assistant",
                    "content": response_text,
                    "usage": usage
                })

            except Exception as e:
                error_msg = f"❌ Error: {str(e)}"
                st.error(error_msg)

                # Log error (audit trail)
                st.session_state.audit_logger.log_error(
                    user_email=user["email"],
                    user_id=user["user_id"],
                    error_type=type(e).__name__,
                    error_message=str(e),
                    servers=st.session_state.selected_servers,
                    session_id=st.session_state.session_id
                )

                st.session_state.messages.append({
                    "role": "assistant",
                    "content": error_msg
                })


def main():
    """Main application."""
    # Require authentication (SSO via OAuth2 Proxy)
    user = require_authentication()

    # Initialize session state
    initialize_session_state()

    # Log session start (first time only)
    if "session_logged" not in st.session_state:
        st.session_state.audit_logger.log_session_start(
            user_email=user["email"],
            user_id=user["user_id"],
            session_id=st.session_state.session_id
        )
        st.session_state.session_logged = True

    # Render sidebar and get settings
    model, max_tokens, show_trace, trace_style = render_sidebar()

    # Main content area
    st.title("🧬 Precision Medicine MCP Chat")
    st.markdown("Chat interface for testing deployed MCP servers on GCP Cloud Run")

    # Check API key
    if not st.session_state.chat_handler:
        st.error("❌ ANTHROPIC_API_KEY not configured")
        st.info("""
        **Setup Instructions:**
        1. Set your API key: `export ANTHROPIC_API_KEY=your_key_here`
        2. Restart this app

        Or create a `.env` file with:
        ```
        ANTHROPIC_API_KEY=your_key_here
        ```
        """)
        return

    # Show server status
    render_server_status()

    st.markdown("---")

    # Render chat history
    render_chat_history(show_trace=show_trace, trace_style=trace_style)

    # Chat input
    # Check for example prompt loaded
    default_value = ""
    if hasattr(st.session_state, 'example_prompt'):
        default_value = st.session_state.example_prompt
        delattr(st.session_state, 'example_prompt')

    if prompt := st.chat_input("Ask me anything about precision medicine...", key="chat_input"):
        handle_user_input(prompt, model, max_tokens)

    # Show example if set
    if default_value:
        handle_user_input(default_value, model, max_tokens)


if __name__ == "__main__":
    main()
