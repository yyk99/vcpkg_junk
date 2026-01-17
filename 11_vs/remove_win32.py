#!/usr/bin/env python3
"""
Script to remove Win32 platform configurations from Visual Studio project files
"""
import re
import sys
from pathlib import Path

def remove_win32_from_vcxproj(file_path):
    """Remove Win32 platform configurations from a .vcxproj file"""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    original_content = content

    # Remove ProjectConfiguration entries for Win32
    content = re.sub(
        r'    <ProjectConfiguration Include="(?:Debug|Release)\|Win32">\s*\n\s*<Configuration>(?:Debug|Release)</Configuration>\s*\n\s*<Platform>Win32</Platform>\s*\n\s*</ProjectConfiguration>\s*\n',
        '',
        content,
        flags=re.MULTILINE
    )

    # Remove PropertyGroup sections for Win32
    content = re.sub(
        r'  <PropertyGroup Condition="[^"]*Win32[^"]*" Label="Configuration">.*?</PropertyGroup>\s*\n',
        '',
        content,
        flags=re.DOTALL
    )

    # Remove ImportGroup sections for Win32
    content = re.sub(
        r'  <ImportGroup Label="PropertySheets" Condition="[^"]*Win32[^"]*">.*?</ImportGroup>\s*\n',
        '',
        content,
        flags=re.DOTALL
    )

    # Remove ItemDefinitionGroup sections for Win32
    content = re.sub(
        r'  <ItemDefinitionGroup Condition="[^"]*Win32[^"]*">.*?</ItemDefinitionGroup>\s*\n',
        '',
        content,
        flags=re.DOTALL
    )

    if content != original_content:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    return False

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python remove_win32.py <vcxproj_file> [<vcxproj_file> ...]")
        sys.exit(1)

    modified_count = 0
    for file_arg in sys.argv[1:]:
        file_path = Path(file_arg)
        if file_path.exists() and file_path.suffix == '.vcxproj':
            if remove_win32_from_vcxproj(file_path):
                print(f"Modified: {file_path}")
                modified_count += 1
            else:
                print(f"No changes needed: {file_path}")
        else:
            print(f"Skipped: {file_path} (not found or not a .vcxproj file)")

    print(f"\nTotal files modified: {modified_count}")
