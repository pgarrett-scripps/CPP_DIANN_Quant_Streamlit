import os
import sys
import subprocess


def run_streamlit_app(module_path):
    """Helper function to run a Streamlit app from a module path."""
    # Get the directory of the current script
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct absolute path if relative
    if not os.path.isabs(module_path):
        module_path = os.path.join(current_dir, module_path)

    # Run streamlit command
    cmd = [sys.executable, "-m", "streamlit", "run", module_path, "--server.maxUploadSize","5000"]
    
    try:
        subprocess.run(cmd)
    except KeyboardInterrupt:
        print("\n\nShutting down gracefully...")
        sys.exit(0)
    except Exception as e:
        print(f"Error running Streamlit app: {e}")
        sys.exit(1)


def run_app():
    try:

        """Run the sage_input_app.py streamlit application."""
        app_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "app.py"
        )
        run_streamlit_app(app_path)
    except KeyboardInterrupt:
        print("\n\nShutting down Sage input app gracefully...")
        sys.exit(0)

