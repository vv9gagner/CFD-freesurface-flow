import type { Metadata } from "next";
import { Geist, Geist_Mono } from "next/font/google";
import Link from "next/link";
import { Suspense } from "react";

import RunForm from "@/components/RunForm";
import RunList from "@/components/RunList";
import { RunsProvider } from "@/components/RunsProvider";
import "./globals.css";

const geistSans = Geist({
  variable: "--font-geist-sans",
  subsets: ["latin"],
});

const geistMono = Geist_Mono({
  variable: "--font-geist-mono",
  subsets: ["latin"],
});

export const metadata: Metadata = {
  title: "Freesurf – free-surface wave solver",
  description: "Run and explore 2D and 3D two-phase Navier–Stokes simulations of sloshing waves.",
};

export default function RootLayout({ children }: LayoutProps<"/">) {
  return (
    <html lang="en" className={`${geistSans.variable} ${geistMono.variable} h-full antialiased`}>
      <body className="min-h-full">
        <RunsProvider>
          <div className="flex min-h-screen flex-col lg:flex-row">
            <aside className="flex w-full shrink-0 flex-col gap-6 border-b border-line bg-surface p-4 lg:sticky lg:top-0 lg:h-screen lg:w-[340px] lg:overflow-y-auto lg:border-r lg:border-b-0">
              <Link href="/" className="block">
                <h1 className="text-lg font-semibold tracking-tight">Freesurf</h1>
                <p className="text-sm text-ink-2">2D &amp; 3D two-phase Navier–Stokes · VOF free surface</p>
              </Link>
              <RunForm />
              <Suspense fallback={<p className="text-sm text-muted">Loading runs…</p>}>
                <RunList />
              </Suspense>
            </aside>
            <main className="min-w-0 flex-1 p-4 lg:p-6">{children}</main>
          </div>
        </RunsProvider>
      </body>
    </html>
  );
}
