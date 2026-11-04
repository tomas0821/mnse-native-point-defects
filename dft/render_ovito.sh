#!/bin/bash
# Run the OVITO renderer headless. OVITO lives in ~/.venvs/ovito; libOpenGL comes from a local unpack of
# libopengl0 (no root), and Vulkan from the system software driver (lavapipe).
L=$(dirname "$(find "$HOME/.venvs/ovito/extra" -name 'libOpenGL.so.0' | head -1)")
exec env -u DISPLAY -u WAYLAND_DISPLAY LD_LIBRARY_PATH="$L" QT_QPA_PLATFORM=ovitoheadless \
  VK_ICD_FILENAMES=/usr/share/vulkan/icd.d/lvp_icd.json \
  "$HOME/.venvs/ovito/bin/python" "$(dirname "$0")/render_ovito_orders.py" "$@"
