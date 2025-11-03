# Development Environment Setup

This guide provides instructions for setting up your development environment for the qlc3d project.

## Prerequisites

- C++ compiler with C++17 support
- CMake
- Git
- Python 3.x (for running tests)

## Setting Up GitHub Copilot in IntelliJ IDEA

GitHub Copilot is an AI-powered code completion tool that can help you write code faster. If you're using IntelliJ IDEA, follow these steps to set it up:

### 1. Install the GitHub Copilot Plugin

1. Open IntelliJ IDEA
2. Go to **File** → **Settings** (on Windows/Linux) or **IntelliJ IDEA** → **Preferences** (on macOS)
3. Select **Plugins** from the left sidebar
4. Click on the **Marketplace** tab
5. Search for "GitHub Copilot"
6. Click **Install** on the GitHub Copilot plugin
7. Restart IntelliJ IDEA when prompted

### 2. Authenticate GitHub Copilot

After installing the plugin, you need to authenticate with your GitHub account:

1. Once IntelliJ IDEA restarts, you should see a notification about GitHub Copilot
2. Click on the notification, or go to **Tools** → **GitHub Copilot** → **Login to GitHub**
3. A dialog will appear with an authentication code (e.g., `XXXX-XXXX`)
4. Click **Copy and Open** or manually copy the code
5. Your browser will open to https://github.com/login/device
6. **Paste the authentication code** in the text field on the GitHub page
7. Click **Continue**
8. Authorize GitHub Copilot to access your GitHub account
9. Return to IntelliJ IDEA - you should now be authenticated

### Alternative Authentication Method

If the automatic flow doesn't work:

1. Go to **File** → **Settings** → **Tools** → **GitHub Copilot**
2. Click **Sign in to GitHub**
3. Follow the same device authentication flow as described above
4. Enter the authentication code at https://github.com/login/device

### Verifying the Setup

To verify that GitHub Copilot is working:

1. Open any C++ file in your project
2. Start typing a function or comment
3. You should see gray text suggestions from Copilot
4. Press **Tab** to accept a suggestion

### Troubleshooting

- If you don't see Copilot suggestions, check that:
  - You have an active GitHub Copilot subscription
  - The Copilot plugin is enabled in IntelliJ IDEA
  - The status bar at the bottom shows the Copilot icon (not disabled)
  
- If authentication fails:
  - Try logging out and logging back in through **Tools** → **GitHub Copilot**
  - Make sure you're using the correct GitHub account with Copilot access
  - Check your internet connection

## Building the Project

For instructions on building qlc3d, see the main [README.md](README.md).

## Running Tests

While in the build directory:
```bash
ctest
```

## Contributing

When contributing to this project, please ensure:
- Your code follows the existing style conventions
- All tests pass before submitting a pull request
- You've tested your changes locally
