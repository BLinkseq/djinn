// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import starlightUiTweaks from 'starlight-ui-tweaks'
import starlightScrollToTop from 'starlight-scroll-to-top'
import starlightSidebarTopics from 'starlight-sidebar-topics'
import codeblocks from 'starlight-codeblocks';

export default defineConfig({
  site: 'https://blinkseq.github.io/djinn',
  base: '/djinn',
	integrations: [
    starlight({
      customCss: ['./src/style/global.css'],
      title: 'Djinn',
      logo: { src: '/src/assets/logo.svg', replacesTitle: true },
      favicon: '/favicon.png',
			social: [{ icon: 'github', label: 'GitHub', href: 'https://github.com/BLinkseq/djinn' }],
			lastUpdated: true,
      pagination: false,
      editLink: {
        baseUrl: 'https://github.com/blinkseq/djinn/edit/main/docs/',
      },
      plugins: [
        codeblocks(),
        starlightUiTweaks({
          navbarLinks: [
           { label : 'Version: 3.0', href: "javascript:void(0);"},
           { label: "BLink-seq", href: "https://blinkseq.github.io"},
        ],
      }),
        starlightScrollToTop({
          showProgressRing: true,
          progressRingColor: '#b288f5',
          borderRadius: '50',
        }),
        starlightSidebarTopics([
          {
            label: 'Install',
            link: '/install',
            icon: 'download',
            items: ['install']
          },
          {
            label: 'FASTQ',
            link: '/fastq/',
            icon: 'information',
            items: [{
              label: 'FASTQ File Operations',
              items: [{ autogenerate: { directory: 'fastq' } }],
              }]
          },
          {
            label: 'SAM',
            link: '/sam/',
            icon: 'information',
            items: [{
              label: 'SAM/BAM File Operations',
              items: [{ autogenerate: { directory: 'sam' } }],
            }],
          }
        ],
        ),
      ],
		}),
	],
});
